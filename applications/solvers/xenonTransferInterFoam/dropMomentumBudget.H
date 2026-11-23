// Per-drop momentum budget (diagnostic): over the cells with dropId == j
// (drop + halo), the momentum change rate and the force of each term of the
// momentum equation, evaluated with the end-of-step fields, per drop mass
// [m/s2]. resid = dPdt - sum of the terms: momentum the discretisation creates.
// Sp = U*(ddt(rho) + div(rhoPhi)) counts only with momentumSp.
// Also the drop velocity by momentum (Uz_mom) and by alpha-weighted U (Uz_alpha).

void dropMomentumBudget
(
    const Foam::volScalarField& alpha,
    const Foam::volScalarField& rho,
    const Foam::volVectorField& U,
    const Foam::surfaceScalarField& rhoPhi,
    const Foam::surfaceScalarField& stf,
    const Foam::surfaceScalarField& ghf,
    const Foam::volScalarField& p_rgh,
    const Foam::volVectorField& visc,   // viscous term as it appears on the LHS
    const bool momentumSp               // Sp term in UEqn: part of the balance
)
{
    using namespace Foam;

    const fvMesh& mesh = U.mesh();
    const auto* idPtr = mesh.findObject<volScalarField>("dropId");
    if (!idPtr) return;
    const scalarField& id = idPtr->primitiveField();

    const label nDrop = returnReduce(label(gMax(id) + 1.5), maxOp<label>());
    if (nDrop <= 0) return;

    const scalarField& V = mesh.V();
    const surfaceScalarField& magSf = mesh.magSf();

    enum { DPDT, CONV, SP, GRAV, PRESS, CSF, VISC, NT };
    const char* names[NT] = {"dPdt", "conv", "Sp", "grav", "press", "csf", "visc"};

    const volScalarField massErr(fvc::ddt(rho) + fvc::div(rhoPhi));

    PtrList<volVectorField> term(NT);
    term.set(DPDT, new volVectorField("dPdt", (rho*U - rho.oldTime()*U.oldTime())/mesh.time().deltaT()));
    term.set(CONV, new volVectorField("conv", -fvc::div(rhoPhi, U)));
    term.set(SP, new volVectorField("Sp", massErr*U));
    term.set(GRAV, new volVectorField("grav", fvc::reconstruct(-ghf*fvc::snGrad(rho)*magSf)));
    term.set(PRESS, new volVectorField("press", fvc::reconstruct(-fvc::snGrad(p_rgh)*magSf)));
    term.set(CSF, new volVectorField("csf", fvc::reconstruct(stf*magSf)));
    term.set(VISC, new volVectorField("visc", -visc));

    List<vectorField> sum(NT, vectorField(nDrop, Zero));
    scalarField mass(nDrop, Zero), vol(nDrop, Zero);
    vectorField mom(nDrop, Zero), aU(nDrop, Zero);

    forAll(id, celli)
    {
        if (id[celli] < -0.5) continue;
        const label j = label(id[celli] + 0.5);
        const scalar v = V[celli];
        for (label k = 0; k < NT; ++k) sum[k][j] += v*term[k][celli];
        mass[j] += v*rho[celli];
        mom[j] += v*rho[celli]*U[celli];
        vol[j] += v*alpha[celli];
        aU[j] += v*alpha[celli]*U[celli];
    }
    for (label k = 0; k < NT; ++k) Pstream::listCombineReduce(sum[k], plusEqOp<vector>());
    Pstream::listCombineReduce(mass, plusEqOp<scalar>());
    Pstream::listCombineReduce(vol, plusEqOp<scalar>());
    Pstream::listCombineReduce(mom, plusEqOp<vector>());
    Pstream::listCombineReduce(aU, plusEqOp<vector>());

    for (label j = 0; j < nDrop; ++j)
    {
        if (mass[j] <= VSMALL || vol[j] <= VSMALL) continue;
        vector rhs(Zero);
        for (label k = CONV; k < NT; ++k)
        {
            if (k != SP || momentumSp) rhs += sum[k][j];
        }

        Info<< "dropMomentum " << j
            << " Uz_mom " << mom[j].z()/mass[j]
            << " Uz_alpha " << aU[j].z()/vol[j]
            << " | z [m/s2]";
        for (label k = 0; k < NT; ++k)
        {
            Info<< ' ' << names[k] << ' ' << sum[k][j].z()/mass[j];
        }
        Info<< " resid " << (sum[DPDT][j].z() - rhs.z())/mass[j]
            << " | resid_xy " << (sum[DPDT][j].x() - rhs.x())/mass[j]
            << ' ' << (sum[DPDT][j].y() - rhs.y())/mass[j] << nl;
    }
}
