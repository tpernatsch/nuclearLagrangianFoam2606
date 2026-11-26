// Per-drop net-force control on the regions dropId == j (drop + halo; the
// dropId field comes from gasDiffusionSimple).
//
// Face forces f are force densities times magSf, as reconstructed into U.
// net(f, w) = sum_c w_c reconstruct(f)_c over drop j. reconstruct is not
// conservative, so a balanced interface (CSF against the Laplace pressure jump)
// leaves a spurious net force. balance() adds -(b_j & (x - xc_j))*snGrad(alpha)
// *magSf, a gravity-form term acting only where snGrad(alpha) does, with b_j
// from the 3x3 system that sets net(f + correction, w) to a target
// (zero, or weight() when f includes the gravity term).

class dropNetForce
{
    const Foam::fvMesh& mesh_;
    Foam::label nDrop_;
    Foam::labelList cid_;                  // cell drop id
    Foam::labelList fid_;                  // internal face drop id (either side)
    Foam::List<Foam::labelList> pown_, pnbr_;  // patch face: owner / neighbour drop id
    Foam::vectorField xc_;
    Foam::scalarField vol_, mass_;
    Foam::scalarField rhoV_, volR_, rhoMin_;   // region (drop + halo) sums, min rho
    Foam::surfaceScalarField snGradAlpha_;

public:

    dropNetForce(const Foam::volScalarField& alpha, const Foam::volScalarField& rho)
    :
        mesh_(alpha.mesh()),
        nDrop_(0),
        snGradAlpha_(Foam::fvc::snGrad(alpha))
    {
        using namespace Foam;

        const auto* idPtr = mesh_.findObject<volScalarField>("dropId");
        if (!idPtr)
        {
            FatalErrorInFunction
                << "per-drop net force control needs the dropId field (gasDiffusionSimple)"
                << exit(FatalError);
        }

        volScalarField idv
        (
            IOobject("dropIdNetForce", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
            *idPtr
        );
        idv.correctBoundaryConditions();

        auto toId = [](const scalar s) { return s < -0.5 ? -1 : label(s + 0.5); };

        nDrop_ = returnReduce(label(gMax(idv.primitiveField()) + 1.5), maxOp<label>());

        cid_.setSize(mesh_.nCells());
        forAll(cid_, celli) cid_[celli] = toId(idv[celli]);

        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        fid_.setSize(mesh_.nInternalFaces());
        forAll(fid_, facei) fid_[facei] = max(cid_[own[facei]], cid_[nei[facei]]);

        pown_.setSize(mesh_.boundary().size());
        pnbr_.setSize(mesh_.boundary().size());
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& pp = mesh_.boundary()[patchi];
            const labelUList& fc = pp.faceCells();
            pown_[patchi].setSize(pp.size());
            pnbr_[patchi].setSize(pp.size(), -1);
            const scalarField nbr
            (
                pp.coupled()
              ? idv.boundaryField()[patchi].patchNeighbourField()()
              : scalarField(pp.size(), -1)
            );
            forAll(fc, i)
            {
                pown_[patchi][i] = cid_[fc[i]];
                pnbr_[patchi][i] = toId(nbr[i]);
            }
        }

        const scalarField& V = mesh_.V();
        const vectorField& C = mesh_.C();
        xc_.setSize(nDrop_, Zero);
        vol_.setSize(nDrop_, Zero);
        mass_.setSize(nDrop_, Zero);
        rhoV_.setSize(nDrop_, Zero);
        volR_.setSize(nDrop_, Zero);
        rhoMin_.setSize(nDrop_, GREAT);
        forAll(cid_, celli)
        {
            const label j = cid_[celli];
            if (j < 0) continue;
            const scalar w = alpha[celli]*V[celli];
            xc_[j] += w*C[celli];
            vol_[j] += w;
            mass_[j] += w*rho[celli];
            rhoV_[j] += rho[celli]*V[celli];
            volR_[j] += V[celli];
            rhoMin_[j] = min(rhoMin_[j], rho[celli]);
        }
        Pstream::listCombineReduce(xc_, plusEqOp<vector>());
        Pstream::listCombineReduce(vol_, plusEqOp<scalar>());
        Pstream::listCombineReduce(mass_, plusEqOp<scalar>());
        Pstream::listCombineReduce(rhoV_, plusEqOp<scalar>());
        Pstream::listCombineReduce(volR_, plusEqOp<scalar>());
        Pstream::listCombineReduce(rhoMin_, minEqOp<scalar>());
        forAll(xc_, j) xc_[j] /= max(vol_[j], VSMALL);
    }

    Foam::label nDrops() const { return nDrop_; }
    const Foam::scalarField& vol() const { return vol_; }
    const Foam::scalarField& mass() const { return mass_; }

    // Region mass (drop + halo, the momentum of the dropMomentumBudget log)
    const Foam::scalarField& regionMass() const { return rhoV_; }

    // Weight minus buoyancy of the region, sum (rho - rho_gas) V g; rho_gas is
    // the region minimum (halo in pure gas)
    Foam::vectorField weight(const Foam::vector& g) const
    {
        Foam::vectorField W(nDrop_, Foam::Zero);
        forAll(W, j)
        {
            if (volR_[j] > Foam::VSMALL) W[j] = (rhoV_[j] - rhoMin_[j]*volR_[j])*g;
        }
        return W;
    }

    Foam::vectorField net(const Foam::surfaceScalarField& f, const Foam::scalarField& w) const
    {
        using namespace Foam;
        const volVectorField Fc(fvc::reconstruct(f));
        vectorField F(nDrop_, Zero);
        forAll(cid_, celli)
        {
            if (cid_[celli] >= 0) F[cid_[celli]] += w[celli]*Fc[celli];
        }
        Pstream::listCombineReduce(F, plusEqOp<vector>());
        return F;
    }

    // -(b_j & (x - xc_j))*snGrad(alpha)*magSf on drop faces, dimensions of f
    Foam::tmp<Foam::surfaceScalarField> correction
    (
        const Foam::vectorField& b,
        const Foam::surfaceScalarField& f
    ) const
    {
        using namespace Foam;
        auto tc = surfaceScalarField::New("dropNetForceCorr", IOobject::NO_REGISTER, 0*f);
        auto& c = tc.ref();
        const vectorField& Cf = mesh_.Cf();
        const scalarField& magSf = mesh_.magSf();
        forAll(fid_, facei)
        {
            const label j = fid_[facei];
            if (j >= 0)
            {
                c[facei] = -(b[j] & (Cf[facei] - xc_[j]))*snGradAlpha_[facei]*magSf[facei];
            }
        }
        auto& cbf = c.boundaryFieldRef();
        forAll(cbf, patchi)
        {
            const vectorField& pCf = mesh_.Cf().boundaryField()[patchi];
            const scalarField& pMagSf = mesh_.magSf().boundaryField()[patchi];
            const scalarField& sg = snGradAlpha_.boundaryField()[patchi];
            forAll(cbf[patchi], i)
            {
                const label j = max(pown_[patchi][i], pnbr_[patchi][i]);
                if (j >= 0) cbf[patchi][i] = -(b[j] & (pCf[i] - xc_[j]))*sg[i]*pMagSf[i];
            }
        }
        return tc;
    }

    // Correction with net(f + correction, w) = target; F0 returns net(f, w)
    Foam::tmp<Foam::surfaceScalarField> balance
    (
        const Foam::surfaceScalarField& f,
        const Foam::scalarField& w,
        const Foam::vectorField& target,
        Foam::vectorField& F0
    ) const
    {
        using namespace Foam;

        F0 = net(f, w);

        List<tensor> M(nDrop_, Zero);
        for (direction k = 0; k < vector::nComponents; ++k)
        {
            vector e(Zero);
            e[k] = 1;
            const vectorField Fk(net(correction(vectorField(nDrop_, e), f), w));
            forAll(M, j)
            {
                for (direction i = 0; i < vector::nComponents; ++i)
                {
                    M[j](i, k) = Fk[j][i];
                }
            }
        }

        vectorField b(nDrop_, Zero);
        forAll(b, j)
        {
            if (vol_[j] > VSMALL && mag(det(M[j])) > VSMALL*pow3(vol_[j]))
            {
                b[j] = inv(M[j]) & (target[j] - F0[j]);
            }
        }
        return correction(b, f);
    }
};
