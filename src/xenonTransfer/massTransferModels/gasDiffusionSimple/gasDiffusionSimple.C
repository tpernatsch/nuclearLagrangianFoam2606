/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Gas-diffusion based interfaceCompositionModel

    - Henry-law equilibrium at the interface
    - Gas-side diffusion with Sherwood correlation
\*---------------------------------------------------------------------------*/

#include "gasDiffusionSimple.H"

#include "constants.H"
#include "cutCellIso.H"
#include "DynamicList.H"
#include "volPointInterpolation.H"
#include "fvcVolumeIntegrate.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

#include "PstreamReduceOps.H"

#include "regionSplit.H"
#include "Pstream.H"

using namespace Foam::constant;

namespace
{
}

// * * * * * * * * * * * *  Constructors  * * * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::gasDiffusionSimple<Thermo, OtherThermo>
::gasDiffusionSimple
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),

    isoAlpha_(dict.getOrDefault<scalar>("isoAlpha", 0.5)),

    alphaBandMin_(dict.getOrDefault<scalar>("alphaBandMin", 1e-6)),

    alphaGasSampleMin_
    (
        dict.getOrDefault<scalar>
        (
            "alphaGasSampleMin",
            1.0 - dict.getOrDefault<scalar>("isoAlpha", 0.5)
        )
    ),

    alphaLiqSampleMin_
    (
        dict.getOrDefault<scalar>
        (
            "alphaLiqSampleMin",
            dict.getOrDefault<scalar>("isoAlpha", 0.5)
        )
    ),

    alphaLo_(dict.getOrDefault<scalar>("alphaLo", 0.05)),
    alphaHi_(dict.getOrDefault<scalar>("alphaHi", 0.95)),
    alphaClusterCut_(dict.getOrDefault<scalar>("alphaClusterCut", 0.5)),
    minDropCells_(dict.getOrDefault<label>("minDropCells", 2)),

    dispersedPhaseName_
    (
        dict.getOrDefault<word>("dispersedPhase", word("liquid"))
    ), 

    geometryVolumeFactor_
    (
        dict.getOrDefault<scalar>("geometryVolumeFactor", 1.0)
    ),

    C_("C", dimless, dict),

    Dgas_
    (
        dict.getOrDefault<scalar>
        (
            "Dgas",
            dict.getOrDefault<scalar>("D_gas", 1e-5)
        )
    ),

    Dliq_
    (
        dict.getOrDefault<scalar>
        (
            "Dliq",
            dict.getOrDefault<scalar>("D_liquid", 2.0e-9)
        )
    ),

    gasSherwoodModel_
    (
        dict.getOrDefault<word>
        (
            "gasSherwoodModel",
            dict.getOrDefault<word>("ShModel", "fixed")
        )
    ),

    gasSherwood_
    (
        dict.getOrDefault<scalar>
        (
            "gasSherwood",
            dict.getOrDefault<scalar>("Sh", 2.0)
        )
    ),

    liquidSherwoodModel_
    (
        dict.getOrDefault<word>
        (
            "liquidSherwoodModel",
            dict.getOrDefault<word>("ShLiquidModel", "fixed")
        )
    ),

    liquidSherwood_
    (
        dict.getOrDefault<scalar>
        (
            "liquidSherwood",
            dict.getOrDefault<scalar>("ShLiquid", 2.0)
        )
    ),

    muGas_(dict.getOrDefault<scalar>("muGas", SMALL)),
    dp_(dict.getOrDefault<scalar>("dp", 0.0)),
    dpMode_(dict.getOrDefault<word>("dpMode", "fixed")),

    // Henry/partition equilibrium:
    k_(dict.get<scalar>("k")),

    debugDeff_(dict.getOrDefault<bool>("debugDeff", false)),
    debugDrops_(dict.getOrDefault<bool>("debugDrops", false)),
    writeDropFields_(dict.getOrDefault<bool>("writeDropFields", false)),

    debugInternalCirculation_
    (
        dict.getOrDefault<bool>("debugInternalCirculation", true)
    ),

    interfaceArea_
    (
        IOobject
        (
            "interfaceArea",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),

    mDotc_
    (
        IOobject
        (
            "mDotc",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity/dimTime, Zero)
    ),

    dropIdField_
    (
        IOobject
        (
            "dropId",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writeDropFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("minusOne", dimless, -1.0)
    ),

    dropDpField_
    (
        IOobject
        (
            "dropDp",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writeDropFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimLength, Zero)
    ),

    dropReField_
    (
        IOobject
        (
            "dropRe",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writeDropFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    dropShGasField_
    (
        IOobject
        (
            "dropShGas",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writeDropFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    dropShLiqField_
    (
        IOobject
        (
            "dropShLiq",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            writeDropFields_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimless, Zero)
    )
{
    isoAlpha_ = max(min(isoAlpha_, scalar(1) - SMALL), SMALL);
    alphaBandMin_ = max(alphaBandMin_, scalar(0));

    alphaGasSampleMin_ =
        min(max(alphaGasSampleMin_, alphaBandMin_), scalar(1));

    alphaLiqSampleMin_ =
        min(max(alphaLiqSampleMin_, alphaBandMin_), scalar(1));

    alphaLo_ = max(min(alphaLo_, scalar(0.49)), scalar(0));
    alphaHi_ = min(max(alphaHi_, scalar(0.51)), scalar(1));
    alphaClusterCut_ = max(min(alphaClusterCut_, scalar(1)), scalar(0));
    minDropCells_ = max(minDropCells_, 1);
    geometryVolumeFactor_ = max(geometryVolumeFactor_, scalar(1.0));

    if (alphaHi_ <= alphaLo_)
    {
        alphaLo_ = 0.05;
        alphaHi_ = 0.95;
    }

    if
    (
        gasSherwoodModel_ != "fixed"
    && gasSherwoodModel_ != "RanzMarshall"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported gasSherwoodModel '" << gasSherwoodModel_
            << "'. Valid choices are fixed | RanzMarshall"
            << exit(FatalIOError);
    }

    if
    (
        liquidSherwoodModel_ != "fixed"
    && liquidSherwoodModel_ != "Higbie"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported liquidSherwoodModel '" << liquidSherwoodModel_
            << "'. Valid choices are fixed | Higbie"
            << exit(FatalIOError);
    }

    if (gasSherwoodModel_ == "RanzMarshall" && muGas_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "gasSherwoodModel RanzMarshall requires positive muGas"
            << exit(FatalIOError);
    }

    if (gasSherwood_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "gasDiffusionSimple requires positive gasSherwood. "
            << "Current gasSherwood = " << gasSherwood_
            << exit(FatalIOError);
    }

    if (liquidSherwood_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "gasDiffusionSimple requires positive liquidSherwood. "
            << "Current liquidSherwood = " << liquidSherwood_
            << exit(FatalIOError);
    }

    if (dpMode_ != "fixed" && dpMode_ != "globalVA" && dpMode_ != "perDropVA")
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported dpMode '" << dpMode_ << "'. Valid choices are fixed | globalVA | perDropVA"
            << exit(FatalIOError);
    }

    if (k_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "gasDiffusionSimple requires k >= 0. "
            << "Current k = " << k_
            << exit(FatalIOError);
    }

    // Dispersed phase accept ("liquid" or "from"), ("gas" or "to")
    if (dispersedPhaseName_ == "from")
    {
        dispersedPhaseName_ = this->pair().from().name();
    }
    else if (dispersedPhaseName_ == "to")
    {
        dispersedPhaseName_ = this->pair().to().name();
    }

    if
    (
        dispersedPhaseName_ != this->pair().from().name()
    && dispersedPhaseName_ != this->pair().to().name()
    )
    {
        FatalIOErrorInFunction(dict)
            << "Invalid dispersedPhase = " << dispersedPhaseName_ << nl
            << "Valid options are: "
            << this->pair().from().name() << ", "
            << this->pair().to().name()
            << ", from, to" << nl
            << exit(FatalIOError);
    }

    if (Dgas_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "gasDiffusionSimple requires positive D_gas. "
            << "Current D_gas = " << Dgas_
            << exit(FatalIOError);
    }

    if (Dliq_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "gasDiffusionSimple requires positive D_liquid. "
            << "Current D_liquid = " << Dliq_
            << exit(FatalIOError);
    }

    if (dpMode_ == "fixed" && dp_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "dpMode fixed requires a positive dp."
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::meltingEvaporationModels::gasDiffusionSimple<Thermo, OtherThermo>
::updateInterface(const volScalarField& T)
{
    (void)T;

    const fvMesh& mesh = this->mesh_;

    const bool dispersedIsFrom =
        (dispersedPhaseName_ == this->pair().from().name());

    const volScalarField& alphaDisp =
        dispersedIsFrom ? this->pair().from() : this->pair().to();

    // from donor cells into receiver cells by DiffusiveMassTransferPhaseSystem.
    const volScalarField& alphaSrc = this->pair().from();

    const volPointInterpolation& vpi = volPointInterpolation::New(mesh);

    scalarField ap(vpi.interpolate(alphaDisp));

    cutCellIso cutCell(mesh, ap);
    const labelListList& cPts = mesh.cellPoints();

    // 1) Iso-based area, stored on the source/donor-side cells
    forAll(interfaceArea_, celli)
    {
        interfaceArea_[celli] = 0.0;

        if (alphaSrc[celli] <= alphaLo_) continue;

        const labelList& pts = cPts[celli];
        if (pts.empty()) continue;

        scalar aMin = GREAT;
        scalar aMax = -GREAT;

        forAll(pts, i)
        {
            const scalar av = ap[pts[i]];
            aMin = min(aMin, av);
            aMax = max(aMax, av);
        }

        const scalar iso  = isoAlpha_;
        const scalar band = max(alphaBandMin_, 0.0);

        if (aMin >= iso - band) continue;
        if (aMax <= iso + band) continue;

        const label status = cutCell.calcSubCell(celli, isoAlpha_);

        if (status == 0)
        {
            interfaceArea_[celli] = mag(cutCell.faceArea())/mesh.V()[celli];
        }
    }

    // 2) Face fallback for sharp interfaces.
    scalarField ifaceFace(mesh.nCells(), 0.0);

    const label nIntFaces = mesh.nInternalFaces();
    const labelList& owner = mesh.owner();
    const labelList& neigh = mesh.neighbour();

    for (label facei = 0; facei < nIntFaces; ++facei)
    {
        const label own = owner[facei];
        const label nei = neigh[facei];

        const scalar aDO = alphaDisp[own];
        const scalar aDN = alphaDisp[nei];

        const bool crossesSplit =
            ((aDO > isoAlpha_) != (aDN > isoAlpha_));

        const bool brackets =
            (min(aDO, aDN) < alphaLo_ && max(aDO, aDN) > alphaHi_);

        if (!(crossesSplit || brackets)) continue;

        const label srcCell =
            (alphaSrc[own] >= alphaSrc[nei]) ? own : nei;

        if (alphaSrc[srcCell] <= alphaLo_) continue;

        const scalar Af = mag(mesh.Sf()[facei]);

        ifaceFace[srcCell] += Af/mesh.V()[srcCell];
    }

    forAll(mesh.boundary(), patchi)
    {
        if (!alphaDisp.boundaryField()[patchi].coupled()) continue;

        const fvPatch& p = mesh.boundary()[patchi];

        const fvPatchScalarField& aDP =
            alphaDisp.boundaryField()[patchi];

        const tmp<scalarField> taDN = aDP.patchNeighbourField();
        const scalarField& aDN = taDN();

        const fvPatchScalarField& aSP =
            alphaSrc.boundaryField()[patchi];

        const tmp<scalarField> taSN = aSP.patchNeighbourField();
        const scalarField& aSN = taSN();

        const labelUList& faceCells = p.faceCells();

        forAll(faceCells, fi)
        {
            const label ownCell = faceCells[fi];

            const scalar aDO = alphaDisp[ownCell];
            const scalar aDNei = aDN[fi];

            const bool crossesSplit =
                ((aDO > isoAlpha_) != (aDNei > isoAlpha_));

            const bool brackets =
                (min(aDO, aDNei) < alphaLo_ && max(aDO, aDNei) > alphaHi_);

            if (!(crossesSplit || brackets)) continue;

            const bool localIsSource = (alphaSrc[ownCell] >= aSN[fi]);

            if (!localIsSource) continue;
            if (alphaSrc[ownCell] <= alphaLo_) continue;

            const scalar Af = mag(p.Sf()[fi]);

            ifaceFace[ownCell] += Af/mesh.V()[ownCell];
        }
    }

    forAll(interfaceArea_, celli)
    {
        if (interfaceArea_[celli] <= SMALL && ifaceFace[celli] > 0.0)
        {
            interfaceArea_[celli] = ifaceFace[celli];
        }
    }
}

// * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::gasDiffusionSimple<Thermo, OtherThermo>::Kexp
(
    const volScalarField& T
)
{
    const fvMesh& mesh = this->mesh_;

    vector gVec(0, 0, -9.81);

    if (mesh.foundObject<uniformDimensionedVectorField>("g"))
    {
        gVec = mesh.lookupObject<uniformDimensionedVectorField>("g").value();
    }

    const scalar gMag = mag(gVec);

    const vector eFall =
        gMag > SMALL ? gVec/gMag : vector(0, 0, -1);

    const scalar volumeFactor = geometryVolumeFactor_;

    // 1) Update interface area field
    updateInterface(T);

    // Species name without phase suffix, e.g. xenon.gas -> xenon
    const word specie(IOobject::member(this->transferSpecie()));

    const Thermo& fromThermo = this->fromThermo_;
    const OtherThermo& toThermo = this->toThermo_;

    const volScalarField& Y_liq = fromThermo.composition().Y(specie);
    const volScalarField& Y_gas = toThermo.composition().Y(specie);

    // Conservative liquid storage field: Zxenon.liquid
    const word zLiqName =
        IOobject::groupName("Z" + specie, this->pair().from().name());

    const volScalarField* Z_liqPtr = nullptr;

    if (mesh.foundObject<volScalarField>(zLiqName))
    {
        Z_liqPtr = &mesh.lookupObject<volScalarField>(zLiqName);
    }
    else if (Pstream::master())
    {
        WarningInFunction
            << "Did not find conservative field " << zLiqName
            << ". Falling back to raw field " << Y_liq.name() << nl;
    }

    // Conservative gas storage field: Zxenon.gas
    const word zGasName =
        IOobject::groupName("Z" + specie, this->pair().to().name());

    const volScalarField* Z_gasPtr = nullptr;

    if (mesh.foundObject<volScalarField>(zGasName))
    {
        Z_gasPtr = &mesh.lookupObject<volScalarField>(zGasName);
    }
    else if (Pstream::master())
    {
        WarningInFunction
            << "Did not find conservative field " << zGasName
            << ". Falling back to raw field " << Y_gas.name() << nl;
    }

    const tmp<volScalarField> trhol = this->pair().from().rho();
    const volScalarField& rhol = trhol();

    const tmp<volScalarField> trhog = this->pair().to().rho();
    const volScalarField& rhog = trhog();

    const tmp<volVectorField> tUliq = this->pair().from().U();
    const volVectorField& Uliq = tUliq();

    const tmp<volVectorField> tUgas = this->pair().to().U();
    const volVectorField& Ugas = tUgas();

    const volScalarField& alpha = this->pair().from();

    // Receiver gas phase.
    const volScalarField& alphaGas = this->pair().to();

    const bool dispersedIsFrom =
        (dispersedPhaseName_ == this->pair().from().name());

    const volScalarField& alphaDisp =
        dispersedIsFrom ? this->pair().from() : this->pair().to();

    const volVectorField& Udisp =
        dispersedIsFrom ? Uliq : Ugas;

    const bool computeInternalCirculation =
        debugInternalCirculation_ && dispersedIsFrom;

    const scalar alphaDivFloor = max(alphaBandMin_, SMALL);
    const scalar alphaGasSampleMin = alphaGasSampleMin_;
    const scalar alphaLiqSampleMin = alphaLiqSampleMin_;

    auto clamp01 = [](const scalar y) -> scalar
    {
        return max(min(y, scalar(1)), scalar(0));
    };

    // 2) Gas-side interfacial sampling
    scalarField YgSum(mesh.nCells(), 0.0);
    scalarField rhogSum(mesh.nCells(), 0.0);
    vectorField UgSum(mesh.nCells(), Zero);
    scalarField wSum(mesh.nCells(), 0.0);

    scalarField alphaGSampleSum(mesh.nCells(), 0.0);
    scalarField YgRawSampleSum(mesh.nCells(), 0.0);
    scalarField ZgSampleSum(mesh.nCells(), 0.0);
    scalarField YgFromZSampleSum(mesh.nCells(), 0.0);

    const label nIntFaces = mesh.nInternalFaces();
    const labelList& owner = mesh.owner();
    const labelList& neigh = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    for (label facei = 0; facei < nIntFaces; ++facei)
    {
        const label own = owner[facei];
        const label nei = neigh[facei];

        const scalar aO = alpha[own];
        const scalar aN = alpha[nei];

        const bool crossesSplit =
            ((aO > isoAlpha_) != (aN > isoAlpha_));

        const bool brackets =
            (min(aO, aN) < alphaLo_ && max(aO, aN) > alphaHi_);

        if (!(crossesSplit || brackets)) continue;

        const label liqCell = (aO >= aN) ? own : nei;
        const label gasCell = (aO >= aN) ? nei : own;

        if (alpha[liqCell] <= alphaLo_) continue;

        const scalar aGSample = alphaGas[gasCell];

        if (aGSample <= alphaGasSampleMin) continue;

        const scalar YgRawSample = clamp01(Y_gas[gasCell]);

        scalar ZgSample = aGSample*YgRawSample;
        scalar YgFromZSample = YgRawSample;
        scalar YgEff = YgRawSample;

        if (Z_gasPtr)
        {
            ZgSample = (*Z_gasPtr)[gasCell];

            YgFromZSample =
                clamp01(ZgSample/max(aGSample, alphaGasSampleMin));

            YgEff = YgFromZSample;
        }

        const scalar w = mag(Sf[facei]);

        YgSum[liqCell] += w*YgEff;
        rhogSum[liqCell] += w*rhog[gasCell];
        UgSum[liqCell] += w*Ugas[gasCell];
        wSum[liqCell] += w;

        alphaGSampleSum[liqCell] += w*aGSample;
        YgRawSampleSum[liqCell] += w*YgRawSample;
        ZgSampleSum[liqCell] += w*ZgSample;
        YgFromZSampleSum[liqCell] += w*YgFromZSample;
    }

    forAll(mesh.boundary(), patchi)
    {
        if (!alpha.boundaryField()[patchi].coupled()) continue;

        const fvPatch& p = mesh.boundary()[patchi];

        const fvPatchScalarField& alphaPatch =
            alpha.boundaryField()[patchi];

        const fvPatchScalarField& YgasPatch =
            Y_gas.boundaryField()[patchi];

        const fvPatchScalarField& rhogPatch =
            rhog.boundaryField()[patchi];

        const fvPatchVectorField& UgasPatch =
            Ugas.boundaryField()[patchi];

        const tmp<scalarField> tAlphaN = alphaPatch.patchNeighbourField();
        const tmp<scalarField> tYgasN = YgasPatch.patchNeighbourField();
        const tmp<scalarField> tRhogN = rhogPatch.patchNeighbourField();
        const tmp<vectorField> tUgasN = UgasPatch.patchNeighbourField();

        const scalarField& alphaN = tAlphaN();
        const scalarField& YgasN = tYgasN();
        const scalarField& rhogN = tRhogN();
        const vectorField& UgasN = tUgasN();

        const fvPatchScalarField& alphaGasPatch =
            alphaGas.boundaryField()[patchi];

        const tmp<scalarField> tAlphaGasN =
            alphaGasPatch.patchNeighbourField();

        const scalarField& alphaGasN = tAlphaGasN();

        tmp<scalarField> tZgasN;
        const scalarField* ZgasNPtr = nullptr;

        if (Z_gasPtr)
        {
            const fvPatchScalarField& ZgasPatch =
                Z_gasPtr->boundaryField()[patchi];

            tZgasN = ZgasPatch.patchNeighbourField();
            ZgasNPtr = &tZgasN();
        }

        const labelUList& faceCells = p.faceCells();
        const vectorField& SfP = mesh.Sf().boundaryField()[patchi];

        forAll(faceCells, fi)
        {
            const label ownCell = faceCells[fi];

            const scalar aO = alpha[ownCell];
            const scalar aN = alphaN[fi];

            const bool crossesSplit =
                ((aO > isoAlpha_) != (aN > isoAlpha_));

            const bool brackets =
                (min(aO, aN) < alphaLo_ && max(aO, aN) > alphaHi_);

            if (!(crossesSplit || brackets)) continue;

            if (aO < aN) continue;
            if (aO <= alphaLo_) continue;

            const scalar aGSample = alphaGasN[fi];

            if (aGSample <= alphaGasSampleMin) continue;

            const scalar YgRawSample = clamp01(YgasN[fi]);

            scalar ZgSample = aGSample*YgRawSample;
            scalar YgFromZSample = YgRawSample;
            scalar YgEff = YgRawSample;

            if (ZgasNPtr)
            {
                ZgSample = (*ZgasNPtr)[fi];

                YgFromZSample =
                    clamp01(ZgSample/max(aGSample, alphaGasSampleMin));

                YgEff = YgFromZSample;
            }

            const scalar w = mag(SfP[fi]);

            YgSum[ownCell] += w*YgEff;
            rhogSum[ownCell] += w*rhogN[fi];
            UgSum[ownCell] += w*UgasN[fi];
            wSum[ownCell] += w;

            alphaGSampleSum[ownCell] += w*aGSample;
            YgRawSampleSum[ownCell] += w*YgRawSample;
            ZgSampleSum[ownCell] += w*ZgSample;
            YgFromZSampleSum[ownCell] += w*YgFromZSample;
        }
    }

    scalarField YgIfc(mesh.nCells(), 0.0);
    scalarField rhogIfc(mesh.nCells(), SMALL);
    vectorField UgIfc(mesh.nCells(), Zero);

    scalarField alphaGSampleIfc(mesh.nCells(), 0.0);
    scalarField YgRawSampleIfc(mesh.nCells(), 0.0);
    scalarField ZgSampleIfc(mesh.nCells(), 0.0);
    scalarField YgFromZSampleIfc(mesh.nCells(), 0.0);

    forAll(wSum, celli)
    {
        if (wSum[celli] > SMALL)
        {
            YgIfc[celli] = YgSum[celli]/wSum[celli];
            rhogIfc[celli] = max(rhogSum[celli]/wSum[celli], SMALL);
            UgIfc[celli] = UgSum[celli]/wSum[celli];

            alphaGSampleIfc[celli] =
                alphaGSampleSum[celli]/wSum[celli];

            YgRawSampleIfc[celli] =
                YgRawSampleSum[celli]/wSum[celli];

            ZgSampleIfc[celli] =
                ZgSampleSum[celli]/wSum[celli];

            YgFromZSampleIfc[celli] =
                YgFromZSampleSum[celli]/wSum[celli];
        }
    }

    // Liquid-side point sampling
    const volPointInterpolation& vpi = volPointInterpolation::New(mesh);

    const scalarField rholP(vpi.interpolate(rhol));
    const scalarField alphaPts(vpi.interpolate(alpha));

    scalarField YliqEffP(vpi.interpolate(Y_liq));

    if (Z_liqPtr)
    {
        const scalarField ZliqP(vpi.interpolate(*Z_liqPtr));

        forAll(YliqEffP, pointi)
        {
            const scalar aLp = alphaPts[pointi];

            if (aLp > alphaLiqSampleMin)
            {
                YliqEffP[pointi] =
                    clamp01(ZliqP[pointi]/max(aLp, alphaDivFloor));
            }
            else
            {
                YliqEffP[pointi] = 0.0;
            }
        }
    }

    const labelListList& cPts = mesh.cellPoints();

    // 3) Parallel-safe droplet clustering
    boolList blockedFace(mesh.nFaces(), false);

    for (label facei = 0; facei < nIntFaces; ++facei)
    {
        const bool ownCore = alphaDisp[owner[facei]] > alphaClusterCut_;
        const bool neiCore = alphaDisp[neigh[facei]] > alphaClusterCut_;

        blockedFace[facei] = (ownCore != neiCore);
    }

    forAll(mesh.boundary(), patchi)
    {
        if (!alphaDisp.boundaryField()[patchi].coupled()) continue;

        const fvPatch& p = mesh.boundary()[patchi];
        const label start = p.start();

        const fvPatchScalarField& aPatch =
            alphaDisp.boundaryField()[patchi];

        const tmp<scalarField> taN = aPatch.patchNeighbourField();
        const scalarField& aN = taN();

        const labelUList& faceCells = p.faceCells();

        forAll(faceCells, fi)
        {
            const bool ownCore =
                alphaDisp[faceCells[fi]] > alphaClusterCut_;

            const bool neiCore =
                aN[fi] > alphaClusterCut_;

            blockedFace[start + fi] = (ownCore != neiCore);
        }
    }

    regionSplit cellRegion(mesh, blockedFace);
    const label nRegions = cellRegion.nRegions();

    labelList regionCoreSize(nRegions, 0);

    forAll(alphaDisp, celli)
    {
        if (alphaDisp[celli] > alphaClusterCut_)
        {
            ++regionCoreSize[cellRegion[celli]];
        }
    }

    forAll(regionCoreSize, regioni)
    {
        reduce(regionCoreSize[regioni], sumOp<label>());
    }

    labelList regionToDrop(nRegions, -1);
    label nDrops = 0;

    for (label regioni = 0; regioni < nRegions; ++regioni)
    {
        if (regionCoreSize[regioni] >= minDropCells_)
        {
            regionToDrop[regioni] = nDrops++;
        }
    }

    labelList dropId(mesh.nCells(), -1);

    forAll(alphaDisp, celli)
    {
        if (alphaDisp[celli] > alphaClusterCut_)
        {
            const label regioni = cellRegion[celli];

            if (regionToDrop[regioni] >= 0)
            {
                dropId[celli] = regionToDrop[regioni];
            }
        }
    }

    const labelListList& cellCells = mesh.cellCells();

    forAll(dropId, celli)
    {
        if (dropId[celli] >= 0) continue;
        if (interfaceArea_[celli] <= SMALL) continue;

        label bestDrop = -1;
        scalar bestAlpha = -GREAT;

        const labelList& nbrs = cellCells[celli];

        forAll(nbrs, ni)
        {
            const label nCell = nbrs[ni];

            if (dropId[nCell] < 0) continue;

            if (alphaDisp[nCell] > bestAlpha)
            {
                bestAlpha = alphaDisp[nCell];
                bestDrop = dropId[nCell];
            }
        }

        if (bestDrop >= 0)
        {
            dropId[celli] = bestDrop;
        }
    }

    // Droplet statistics
    scalarList dropV(nDrops, 0.0);        // old assigned volume, for Dvol comparison
    scalarList dropVarea(nDrops, 0.0);    // area-consistent volume, for dp = 6V/A
    scalarList dropVcore(nDrops, 0.0);    // optional diagnostic
    scalarList dropA(nDrops, 0.0);

    scalarList dropAreaWeight(nDrops, 0.0);
    scalarList dropRhogInt(nDrops, 0.0);
    List<vector> dropUgInt(nDrops, Zero);
    List<vector> dropUlInt(nDrops, Zero);

    List<vector> dropXInt(nDrops, Zero);
    List<vector> dropUcmInt(nDrops, Zero);

    List<vector> dropC(nDrops, Zero);
    List<vector> dropUcm(nDrops, Zero);
    

    List<symmTensor> dropCovInt(nDrops, Zero);

    scalarList dropDfall(nDrops, 0.0);
    scalarList dropDperp(nDrops, 0.0);
    scalarList dropAR(nDrops, 1.0);
    scalarList dropARSigned(nDrops, 1.0);
    scalarList dropDef(nDrops, 0.0);
    scalarList dropDvol(nDrops, 0.0);       // from old dropV
    scalarList dropDvolArea(nDrops, 0.0);   // from Varea, diagnostic only
    scalarList dropDvolCore(nDrops, 0.0);   // from Vcore, diagnostic only

    scalarList dropUrel(nDrops, 0.0);

    scalarList dropDpOld(nDrops, 0.0);      // old dp = 6*dropV/dropA
    scalarList dropDpArea(nDrops, 0.0);     // new dp = 6*dropVarea/dropA
    scalarList dropDp(nDrops, 0.0);         // active dp used by model

    scalarList dropRe(nDrops, 0.0);

    scalarList dropShGas(nDrops, gasSherwood_);
    scalarList dropShLiq(nDrops, liquidSherwood_);

    scalarList dropKg(nDrops, 0.0);
    scalarList dropKLiq(nDrops, 0.0);
    scalarList dropKoverall(nDrops, 0.0);

    scalarList dropFlux(nDrops, 0.0);

    //     apparent liquid-side Sherwood number based on the bulk droplet
    scalarList dropBulkWeight(nDrops, 0.0);       // [kg] liquid mass in owned drop
    scalarList dropYbulkInt(nDrops, 0.0);         // [kg] * [-]
    scalarList dropMxeBulk(nDrops, 0.0);          // [kg] xenon inventory in drop

    scalarList dropAreaEff(nDrops, 0.0);          // [m2]
    scalarList dropYliqAreaInt(nDrops, 0.0);      // [m2] * [-]
    scalarList dropYeqAreaInt(nDrops, 0.0);       // [m2] * [-]
    scalarList dropRholAreaInt(nDrops, 0.0);      // [m2] * [kg/m3]
    scalarList dropDriveAreaInt(nDrops, 0.0);     // [m2] * [-]

    scalarList dropYbulk(nDrops, 0.0);
    scalarList dropYeqAvg(nDrops, 0.0);
    scalarList dropYliqAreaAvg(nDrops, 0.0);
    scalarList dropDriveBulk(nDrops, 0.0);
    scalarList dropDriveAreaAvg(nDrops, 0.0);

    scalarList dropKliqApp(nDrops, 0.0);          // [m/s]
    scalarList dropShLiqApp(nDrops, 0.0);         // [-]

    scalarList dropBetaExcess(nDrops, 0.0);       // [1/s]
    scalarList dropTauExcess(nDrops, GREAT);      // [s]

    scalarList dropFracRateInventory(nDrops, 0.0);// [1/s]
    scalarList dropTauInventory(nDrops, GREAT);   // [s]

    scalarList dropPeLiq(nDrops, 0.0);            // [-]

    scalarList dropUcirc2Int(nDrops, 0.0);
    scalarList dropUcircRms(nDrops, 0.0);
    scalarList dropUcircOverUfall(nDrops, 0.0);
    scalarList dropUcircOverUrel(nDrops, 0.0);
    scalarList dropPeInt(nDrops, 0.0);
    scalarList dropUcirc2Sym(nDrops, 0.0);
    scalarList dropUcircSymOverUrel(nDrops, 0.0);
    scalarList dropUcircRmsSym(nDrops, 0.0);
    scalarList dropPeIntSym(nDrops, 0.0);

    forAll(dropId, celli)
    {
        const label j = dropId[celli];
        if (j < 0) continue;

        const scalar Vcell = mesh.V()[celli];
        const scalar wVol = alphaDisp[celli]*Vcell;
        const scalar localA = interfaceArea_[celli]*Vcell;

        const bool isCoreCell = alphaDisp[celli] > alphaClusterCut_;
        const bool isAreaCell = interfaceArea_[celli] > SMALL;

        if (dispersedIsFrom)
        {
            const scalar aL = max(alpha[celli], scalar(0));

            if (aL > alphaLiqSampleMin)
            {
                scalar YbulkCell = clamp01(Y_liq[celli]);

                if (Z_liqPtr)
                {
                    YbulkCell =
                        clamp01((*Z_liqPtr)[celli]/max(aL, alphaDivFloor));
                }

                const scalar mL =
                    max(rhol[celli], SMALL)*aL*Vcell;

                dropBulkWeight[j] += mL;
                dropYbulkInt[j] += mL*YbulkCell;
                dropMxeBulk[j] += mL*YbulkCell;
            }
        }

        dropV[j] += wVol;

        if (isCoreCell || isAreaCell)
        {
            dropVarea[j] += wVol;
        }

        if (isCoreCell)
        {
            dropVcore[j] += wVol;
        }

        dropA[j] += localA;

        dropXInt[j] += wVol*mesh.C()[celli];
        dropUcmInt[j] += wVol*Udisp[celli];

        if (isAreaCell && wSum[celli] > SMALL)
        {
            dropAreaWeight[j] += localA;
            dropRhogInt[j] += localA*rhogIfc[celli];
            dropUgInt[j] += localA*UgIfc[celli];
            dropUlInt[j] += localA*Uliq[celli];
        }
    }

    for (label j = 0; j < nDrops; ++j)
    {
        reduce(dropV[j], sumOp<scalar>());
        reduce(dropVarea[j], sumOp<scalar>());
        reduce(dropVcore[j], sumOp<scalar>());
        reduce(dropA[j], sumOp<scalar>());

        reduce(dropAreaWeight[j], sumOp<scalar>());
        reduce(dropRhogInt[j], sumOp<scalar>());
        reduce(dropUgInt[j], sumOp<vector>());
        reduce(dropUlInt[j], sumOp<vector>());

        reduce(dropXInt[j], sumOp<vector>());
        reduce(dropUcmInt[j], sumOp<vector>());

        reduce(dropMxeBulk[j], sumOp<scalar>());

        reduce(dropYbulkInt[j], sumOp<scalar>());
        reduce(dropBulkWeight[j], sumOp<scalar>());
    }

    for (label j = 0; j < nDrops; ++j)
    {
        if (dropV[j] > SMALL)
        {
            dropC[j] = dropXInt[j]/dropV[j];
            dropUcm[j] = dropUcmInt[j]/dropV[j];

            const scalar Vphys = volumeFactor*dropV[j];

            dropDvol[j] =
                pow
                (
                    max(6.0*Vphys/mathematical::pi, SMALL),
                    scalar(1.0/3.0)
                );
        }

        if (dropVarea[j] > SMALL)
        {
            const scalar VareaPhys = volumeFactor*dropVarea[j];

            dropDvolArea[j] =
                pow
                (
                    max(6.0*VareaPhys/mathematical::pi, SMALL),
                    scalar(1.0/3.0)
                );
        }

        if (dropVcore[j] > SMALL)
        {
            const scalar VcorePhys = volumeFactor*dropVcore[j];

            dropDvolCore[j] =
                pow
                (
                    max(6.0*VcorePhys/mathematical::pi, SMALL),
                    scalar(1.0/3.0)
                );
        }
    }

    if (computeInternalCirculation)
    {
        forAll(dropId, celli)
        {
            const label j = dropId[celli];
            if (j < 0) continue;

            const scalar wVol = alphaDisp[celli]*mesh.V()[celli];
            if (wVol <= SMALL) continue;

            const vector uPrime = Udisp[celli] - dropUcm[j];
            dropUcirc2Int[j] += wVol*magSqr(uPrime);

            const vector Utrans = (dropUcm[j] & eFall)*eFall;
            const vector uPrimeSym = Udisp[celli] - Utrans;
            dropUcirc2Sym[j] += wVol*magSqr(uPrimeSym);
        }

        for (label j = 0; j < nDrops; ++j)
        {
            reduce(dropUcirc2Int[j], sumOp<scalar>());
            reduce(dropUcirc2Sym[j], sumOp<scalar>());
        }

        for (label j = 0; j < nDrops; ++j)
        {
            if (dropV[j] > SMALL)
            {
                dropUcircRms[j] =
                    sqrt(max(dropUcirc2Int[j]/dropV[j], scalar(0)));

                dropUcircRmsSym[j] =
                    sqrt(max(dropUcirc2Sym[j]/dropV[j], scalar(0)));

                const scalar UfallMag = mag(dropUcm[j] & eFall);

                dropUcircOverUfall[j] =
                    dropUcircRms[j]/max(UfallMag, SMALL);
            }
        }
    }

    forAll(dropId, celli)
    {
        const label j = dropId[celli];
        if (j < 0) continue;
        if (dropV[j] <= SMALL) continue;

        const scalar wVol = alphaDisp[celli]*mesh.V()[celli];
        const vector dx = mesh.C()[celli] - dropC[j];

        dropCovInt[j] += wVol*symmTensor
        (
            dx.x()*dx.x(),
            dx.x()*dx.y(),
            dx.x()*dx.z(),
            dx.y()*dx.y(),
            dx.y()*dx.z(),
            dx.z()*dx.z()
        );
    }

    for (label j = 0; j < nDrops; ++j)
    {
        reduce(dropCovInt[j], sumOp<symmTensor>());
    }

    for (label j = 0; j < nDrops; ++j)
    {
        if (dropV[j] <= SMALL) continue;

        const symmTensor Ccov = dropCovInt[j]/dropV[j];

        const vector Ce = Ccov & eFall;

        const scalar varFall =
            max(eFall & Ce, scalar(0));

        scalar aFall = sqrt(max(5.0*varFall, scalar(0)));

        if (aFall <= SMALL)
        {
            aFall = max(0.5*dropDvol[j], SMALL);
        }

        const scalar Vphys = volumeFactor*dropV[j];

        const scalar aPerp =
            sqrt
            (
                max
                (
                    3.0*Vphys
                   /(4.0*mathematical::pi*max(aFall, SMALL)),
                    scalar(0)
                )
            );

        dropDfall[j] = 2.0*aFall;
        dropDperp[j] = 2.0*aPerp;

        dropARSigned[j] = aFall/max(aPerp, SMALL);

        dropAR[j] =
            max(aFall, aPerp)/max(min(aFall, aPerp), SMALL);

        dropDef[j] =
            (aFall - aPerp)/max(aFall + aPerp, SMALL);
    }

    for (label j = 0; j < nDrops; ++j)
    {
        if (dropA[j] > SMALL)
        {
            dropDpOld[j] = max(6.0*dropV[j]/dropA[j], SMALL);

            dropDpArea[j] = max(6.0*dropVarea[j]/dropA[j], SMALL);

            dropDp[j] = dropDpArea[j];
        }
        else if (dropDvol[j] > SMALL)
        {
            dropDpOld[j] = dropDvol[j];
            dropDpArea[j] = dropDvol[j];
            dropDp[j] = dropDvol[j];
        }

        if (dropAreaWeight[j] > SMALL)
        {
            const scalar rhogAvg =
                max(dropRhogInt[j]/dropAreaWeight[j], SMALL);

            const vector UgAvg = dropUgInt[j]/dropAreaWeight[j];
            const vector UlAvg = dropUlInt[j]/dropAreaWeight[j];

            const scalar Urel = mag(UgAvg - UlAvg);

            dropUrel[j] = Urel;

            const scalar dpEval =
                (dpMode_ == "fixed" && dp_ > SMALL)
            ? max(dp_, SMALL)
            : max(dropDp[j], SMALL);

            if (gasSherwoodModel_ == "RanzMarshall")
            {
                const scalar ReGas =
                    max(rhogAvg*Urel*dpEval/muGas_, scalar(0));

                const scalar ScGas =
                    max(muGas_/(rhogAvg*max(Dgas_, SMALL)), SMALL);

                dropRe[j] = ReGas;

                dropShGas[j] =
                    2.0 + 0.6*sqrt(ReGas)*pow(ScGas, scalar(1.0/3.0));
            }
            else
            {
                dropRe[j] =
                    (muGas_ > SMALL)
                ? rhogAvg*Urel*dpEval/muGas_
                : 0.0;

                dropShGas[j] = gasSherwood_;
            }

            if (liquidSherwoodModel_ == "Higbie")
            {
                const scalar PeLiq =
                    max(Urel*dpEval/max(Dliq_, SMALL), scalar(0));

                dropShLiq[j] =
                    max
                    (
                        liquidSherwood_,   // lower bound, usually 2
                        2.0*sqrt(PeLiq/mathematical::pi)
                    );
            }
            else
            {
                dropShLiq[j] = liquidSherwood_;
            }
        }
        else
        {
            dropShGas[j] = gasSherwood_;
            dropShLiq[j] = liquidSherwood_;
        }
    }

    if (computeInternalCirculation)
    {
        for (label j = 0; j < nDrops; ++j)
        {
            if (dropDp[j] > SMALL)
            {
                dropPeInt[j] =
                    dropUcircRms[j]*dropDp[j]/max(Dliq_, SMALL);

                dropPeIntSym[j] =
                    dropUcircRmsSym[j]*dropDp[j]/max(Dliq_, SMALL);
            }

            if (dropUrel[j] > SMALL)
            {
                dropUcircOverUrel[j] =
                    dropUcircRms[j]/max(dropUrel[j], SMALL);

                dropUcircSymOverUrel[j] =
                    dropUcircRmsSym[j]/max(dropUrel[j], SMALL);
            }
        }
    }

    for (label j = 0; j < nDrops; ++j)
    {
        scalar dpEval = dropDp[j];

        if (dpMode_ == "fixed" && dp_ > SMALL)
        {
            dpEval = dp_;
        }

        if (dpEval > SMALL)
        {
            dropKg[j] =
                max(dropShGas[j], SMALL)*Dgas_/max(dpEval, SMALL);

            dropKLiq[j] =
                max(dropShLiq[j], SMALL)*Dliq_/max(dpEval, SMALL);

            dropKoverall[j] =
                1.0/
                (
                    1.0/max(dropKLiq[j], SMALL)
                + k_/max(dropKg[j], SMALL)
                );
        }
    }

    // 4) Global fallback geometry + assemble mDotc
    scalar dpGeomGlobal = max(dp_, SMALL);

    if (dp_ <= SMALL || dpMode_ == "globalVA" || dpMode_ == "perDropVA")
    {
        scalar VareaGlobalLocal = 0.0;
        scalar AtotIsoLocal = 0.0;

        forAll(alphaDisp, celli)
        {
            const scalar Vcell = mesh.V()[celli];

            const bool isCoreCell = alphaDisp[celli] > alphaClusterCut_;
            const bool isAreaCell = interfaceArea_[celli] > SMALL;

            if (isCoreCell || isAreaCell)
            {
                VareaGlobalLocal += alphaDisp[celli]*Vcell;
            }

            AtotIsoLocal += interfaceArea_[celli]*Vcell;
        }

        scalar VareaGlobal = VareaGlobalLocal;
        scalar AtotIsoMasked = AtotIsoLocal;

        reduce(VareaGlobal, sumOp<scalar>());
        reduce(AtotIsoMasked, sumOp<scalar>());

        if (VareaGlobal > SMALL && AtotIsoMasked > SMALL)
        {
            dpGeomGlobal = max(6.0*VareaGlobal/AtotIsoMasked, SMALL);
        }
    }

    scalar ShGasGlobal = gasSherwood_;
    scalar ShLiqGlobal = liquidSherwood_;

    if (gasSherwoodModel_ == "RanzMarshall" || liquidSherwoodModel_ == "Higbie")
    {
        scalar Atot = 0.0;
        scalar rhogInt = 0.0;
        vector UgInt(Zero);
        vector UlInt(Zero);

        forAll(interfaceArea_, celli)
        {
            if (interfaceArea_[celli] <= SMALL || wSum[celli] <= SMALL)
            {
                continue;
            }

            const scalar localA = interfaceArea_[celli]*mesh.V()[celli];

            Atot += localA;
            rhogInt += localA*rhogIfc[celli];
            UgInt += localA*UgIfc[celli];
            UlInt += localA*Uliq[celli];
        }

        reduce(Atot, sumOp<scalar>());
        reduce(rhogInt, sumOp<scalar>());
        reduce(UgInt, sumOp<vector>());
        reduce(UlInt, sumOp<vector>());

        if (Atot > SMALL)
        {
            const scalar rhogAvg = max(rhogInt/Atot, SMALL);
            const scalar Urel = mag(UgInt/Atot - UlInt/Atot);

            if (gasSherwoodModel_ == "RanzMarshall")
            {
                const scalar ReGas =
                    max(rhogAvg*Urel*dpGeomGlobal/muGas_, scalar(0));

                const scalar ScGas =
                    max(muGas_/(rhogAvg*max(Dgas_, SMALL)), SMALL);

                ShGasGlobal =
                    2.0 + 0.6*sqrt(ReGas)*pow(ScGas, scalar(1.0/3.0));
            }

            if (liquidSherwoodModel_ == "Higbie")
            {
                const scalar PeLiq =
                    max(Urel*dpGeomGlobal/max(Dliq_, SMALL), scalar(0));

                ShLiqGlobal =
                    max
                    (
                        liquidSherwood_,
                        2.0*sqrt(PeLiq/mathematical::pi)
                    );
            }
        }
    }

    const scalar kgGlobal =
        ShGasGlobal*Dgas_/max(dpGeomGlobal, SMALL);

    const scalar kLiqGlobal =
        ShLiqGlobal*Dliq_/max(dpGeomGlobal, SMALL);

    const scalar KoverallGlobal =
        1.0/
        (
            1.0/max(kLiqGlobal, SMALL)
        + k_/max(kgGlobal, SMALL)
        );

    mDotc_ =
        dimensionedScalar("zeroMdot", dimDensity/dimTime, 0.0);

    dropIdField_ =
        dimensionedScalar("minusOne", dimless, -1.0);

    dropDpField_ =
        dimensionedScalar("zeroDp", dimLength, 0.0);

    dropReField_ =
        dimensionedScalar("zero", dimless, 0.0);

    dropShGasField_ =
    dimensionedScalar("zero", dimless, 0.0);

    dropShLiqField_ =
        dimensionedScalar("zero", dimless, 0.0);

    forAll(dropId, celli)
    {
        const label j = dropId[celli];

        if (j < 0) continue;

        dropIdField_[celli] = scalar(j);
        dropDpField_[celli] = dropDp[j];
        dropReField_[celli] = dropRe[j];

        dropShGasField_[celli] = dropShGas[j];
        dropShLiqField_[celli] = dropShLiq[j];
    }

    forAll(mDotc_, celli)
    {
        const scalar A = interfaceArea_[celli];

        if (A <= SMALL) continue;
        if (wSum[celli] <= SMALL) continue;

        const labelList& pts = cPts[celli];

        if (pts.empty()) continue;

        scalar YliqCell = 0.0;
        scalar rholCell = 0.0;
        scalar wLiqPts = 0.0;

        forAll(pts, pi)
        {
            const label p = pts[pi];

            if (alphaPts[p] > alphaLiqSampleMin)
            {
                YliqCell += YliqEffP[p];
                rholCell += rholP[p];
                wLiqPts += 1.0;
            }
        }

        if (wLiqPts <= SMALL) continue;

        YliqCell /= wLiqPts;
        rholCell /= wLiqPts;

        YliqCell = clamp01(YliqCell);
        rholCell = max(rholCell, SMALL);

        const scalar YgCell = YgIfc[celli];
        const scalar rhogCell = rhogIfc[celli];

        scalar YeqCell = k_*(rhogCell/rholCell)*YgCell;

        YeqCell = clamp01(YeqCell);

        const label j = dropId[celli];

        scalar dpCell = dpGeomGlobal;

        scalar ShGasCell = ShGasGlobal;
        scalar ShLiqCell = ShLiqGlobal;

        if (dpMode_ == "fixed" && dp_ > SMALL)
        {
            dpCell = max(dp_, SMALL);
        }
        else if (dpMode_ == "globalVA")
        {
            dpCell = dpGeomGlobal;
        }
        else if (dpMode_ == "perDropVA")
        {

            if (j >= 0 && dropDp[j] > SMALL)
            {
                dpCell = dropDp[j];
                ShGasCell = dropShGas[j];
                ShLiqCell = dropShLiq[j];
            }
        }

        // Gas-side coefficient from gas-side Ranz-Marshall
        const scalar kgCell =
            ShGasCell*Dgas_/max(dpCell, SMALL);

        const scalar kLiqCell =
            ShLiqCell*Dliq_/max(dpCell, SMALL);

        const scalar KoverallCell =
            1.0/
            (
                1.0/max(kLiqCell, SMALL)
            + k_/max(kgCell, SMALL)
            );

        mDotc_[celli] =
            C_.value()
        *KoverallCell
        *A
        *rholCell
        *(YliqCell - YeqCell);

        if (j >= 0)
        {
            const scalar cellFlux =
                mDotc_[celli]*mesh.V()[celli];

            dropFlux[j] += cellFlux;

            const scalar localAphys =
                A*mesh.V()[celli];

            if (localAphys > SMALL)
            {
                dropAreaEff[j] += localAphys;
                dropYliqAreaInt[j] += localAphys*YliqCell;
                dropYeqAreaInt[j] += localAphys*YeqCell;
                dropRholAreaInt[j] += localAphys*rholCell;
                dropDriveAreaInt[j] += localAphys*(YliqCell - YeqCell);
            }
        }
    }

    for (label j = 0; j < nDrops; ++j)
    {
        reduce(dropFlux[j], sumOp<scalar>());

        reduce(dropAreaEff[j], sumOp<scalar>());
        reduce(dropYliqAreaInt[j], sumOp<scalar>());
        reduce(dropYeqAreaInt[j], sumOp<scalar>());
        reduce(dropRholAreaInt[j], sumOp<scalar>());
        reduce(dropDriveAreaInt[j], sumOp<scalar>());
    }

    for (label j = 0; j < nDrops; ++j)
    {
        if (dropDp[j] > SMALL)
        {
            dropPeLiq[j] =
                dropUrel[j]*dropDp[j]/max(Dliq_, SMALL);
        }

        if
        (
            dispersedIsFrom
        && dropBulkWeight[j] > SMALL
        && dropAreaEff[j] > SMALL
        && dropDp[j] > SMALL
        )
        {
            dropYbulk[j] =
                dropYbulkInt[j]/dropBulkWeight[j];

            dropYliqAreaAvg[j] =
                dropYliqAreaInt[j]/dropAreaEff[j];

            dropYeqAvg[j] =
                dropYeqAreaInt[j]/dropAreaEff[j];

            const scalar rholAreaAvg =
                max(dropRholAreaInt[j]/dropAreaEff[j], SMALL);

            dropDriveBulk[j] =
                dropYbulk[j] - dropYeqAvg[j];

            dropDriveAreaAvg[j] =
                dropDriveAreaInt[j]/dropAreaEff[j];

            if (mag(dropDriveBulk[j]) > SMALL)
            {
                dropKliqApp[j] =
                    dropFlux[j]
                /(
                        max(dropAreaEff[j], SMALL)
                    *rholAreaAvg
                    *dropDriveBulk[j]
                    );

                dropShLiqApp[j] =
                    dropKliqApp[j]*dropDp[j]/max(Dliq_, SMALL);

                dropBetaExcess[j] =
                    dropFlux[j]
                /(
                        max(dropBulkWeight[j], SMALL)
                    *dropDriveBulk[j]
                    );

                dropTauExcess[j] =
                    scalar(1)/max(mag(dropBetaExcess[j]), SMALL);
            }

            if (dropMxeBulk[j] > SMALL)
            {
                dropFracRateInventory[j] =
                    dropFlux[j]/dropMxeBulk[j];

                dropTauInventory[j] =
                    scalar(1)/max(mag(dropFracRateInventory[j]), SMALL);
            }
        }
    }

    scalar mDotIntegral = 0.0;

    forAll(mDotc_, celli)
    {
        mDotIntegral += mDotc_[celli]*mesh.V()[celli];
    }

    reduce(mDotIntegral, sumOp<scalar>());

    const scalar mDotMin = gMin(mDotc_);
    const scalar mDotMax = gMax(mDotc_);

    // 5) Compact diagnostics
    if ((debugDeff_ || debugDrops_) && Pstream::master())
    {
        Info<< "gasDiffusionSimple:"
            << " dispersedPhase=" << dispersedPhaseName_
            << " nDrops=" << nDrops
            << " dpMode=" << dpMode_
            << " dpGlobal=" << dpGeomGlobal
            << " mDotcIntegral=" << mDotIntegral
            << " mDotcMin=" << mDotMin
            << " mDotcMax=" << mDotMax
            << " Dgas=" << Dgas_
            << " Dliq=" << Dliq_
            << " gasSherwoodModel=" << gasSherwoodModel_
            << " ShGasGlobal=" << ShGasGlobal
            << " liquidSherwoodModel=" << liquidSherwoodModel_
            << " ShLiqGlobal=" << ShLiqGlobal
            << " kgGlobal=" << kgGlobal
            << " kLiqGlobal=" << kLiqGlobal
            << " KoverallGlobal=" << KoverallGlobal
            << nl;
    }

    if (debugDrops_ && Pstream::master())
    {
        for (label j = 0; j < nDrops; ++j)
        {
            const scalar Aphys =
                volumeFactor*dropA[j];

            const scalar Asphere =
                mathematical::pi*sqr(dropDvol[j]);

            const scalar areaRatio =
                Aphys/max(Asphere, SMALL);

            const scalar DvolOverDp =
                dropDvol[j]/max(dropDp[j], SMALL);

            const scalar VcoreOverV =
                dropVcore[j]/max(dropV[j], SMALL);

            Info<< "  drop " << j
                << ": V=" << dropV[j]
                << " Vcore=" << dropVcore[j]
                << " A=" << dropA[j]
                << " dp=" << dropDp[j]
                << " Dvol=" << dropDvol[j]
                << " Dfall=" << dropDfall[j]
                << " Dperp=" << dropDperp[j]
                << " areaRatio=" << areaRatio
                << " DvolOverDp=" << DvolOverDp
                << " VcoreOverV=" << VcoreOverV
                << " AR=" << dropAR[j]
                << " ARsigned=" << dropARSigned[j]
                << " deformation=" << dropDef[j]
                << " Ccm=" << dropC[j]
                << " Ucm=" << dropUcm[j]
                << " Ufall=" << (dropUcm[j] & eFall)
                << " Urel=" << dropUrel[j]
                << " Re=" << dropRe[j];

            if (computeInternalCirculation)
            {
                Info<< " UcircRms=" << dropUcircRms[j]
                    << " UcircOverUfall=" << dropUcircOverUfall[j]
                    << " UcircOverUrel=" << dropUcircOverUrel[j]
                    << " PeInt=" << dropPeInt[j]
                    << " UcircRmsSym=" << dropUcircRmsSym[j]
                    << " UcircSymOverUrel=" << dropUcircSymOverUrel[j]
                    << " PeIntSym=" << dropPeIntSym[j];
            }

            Info<< " ShGas=" << dropShGas[j]
                << " ShLiq=" << dropShLiq[j]
                << " kg=" << dropKg[j]
                << " kLiqModel=" << dropKLiq[j]
                << " Koverall=" << dropKoverall[j]
                << " integratedFlux=" << dropFlux[j]
                << " MxeDrop=" << dropMxeBulk[j]
                << " Ybulk=" << dropYbulk[j]
                << " YliqAreaAvg=" << dropYliqAreaAvg[j]
                << " YeqAvg=" << dropYeqAvg[j]
                << " driveBulk=" << dropDriveBulk[j]
                << " driveAreaAvg=" << dropDriveAreaAvg[j]
                << " betaExcess=" << dropBetaExcess[j]
                << " tauExcess=" << dropTauExcess[j]
                << " fracRateInventory=" << dropFracRateInventory[j]
                << " tauInventory=" << dropTauInventory[j]
                << " PeLiq=" << dropPeLiq[j]
                << " kLiqApp=" << dropKliqApp[j]
                << " ShLiqApp=" << dropShLiqApp[j]
                << nl;
        }
    }

    return tmp<volScalarField>::New(mDotc_);
}

