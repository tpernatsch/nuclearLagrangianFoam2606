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

/* Example phaseProperties entry

(liquid to gas)
{
    type                    gasDiffusionSimple;
    species                 xenon.gas;

    C                       1;
    dispersedPhase          gas;
    geometryVolumeFactor    4;

    isoAlpha                0.5;
    alphaBandMin            1e-6;
    alphaLo                 0.05;
    alphaHi                 0.95;
    alphaClusterCut         0.5;
    minDropCells            4;
    dropHaloLayers          0;             // attach smeared dispersed phase within N face layers to its drop

    debugDeff               true;
    debugDrops              true;
    debugInternalCirculation false;
    writeDropFields         false;

    k                       2.08e-4;
    Dgas                    1e-5;
    Dliq                    2.006e-9;
    muGas                   4.6913e-05;

    gasSherwoodModel        RanzMarshall;
    gasSherwood             2;
    liquidSherwoodModel     Higbie;        // fixed | Higbie | Newman
    liquidSherwood          2;             // fixed value / lower bound
    liquidDrive             interface;     // interface | bulk (default bulk for Newman)

    dpMode                  perDropVA;     // fixed | globalVA | perDropVA | perDropVol
    dp                      3.0e-4;

    dropArea                iso;           // iso | volumeSphere
    correctDropArea         false;         // spread each drop's area over its flux cells
    slipVelocity            interface;     // interface | shell
    slipShell               (1.5 3);       // shell radii / drop radius (slipVelocity shell)

    contactAge              false;         // transport contactAge.<dispersedPhase> [s]
    ageAlphaCut             1e-3;          // no age below this dispersed-phase fraction
    surfaceRenewal          none;          // none | stretching
}

Contact age theta: time the dispersed phase has spent in drops (clusters
with interface), transported conservatively as alpha*theta with the phase
flux, so it follows drops through series, breakup and mapped re-entry:
    d(alpha*theta)/dt + div(alphaPhi*theta) = alpha*chi_drop - alpha*theta*s
s (surfaceRenewal stretching) is the per-drop interface stretching rate
max(0, int(div_s U) dA / A): new interface is fresh, so the mean surface age
drops. Inflow takes the patch value of theta (0 unless set, e.g. mapped).

liquidSherwoodModel Newman (needs contactAge, liquid drops in gas): stagnant
drop with zero surface concentration (Newman 1931, Crank eq. 6.20). The
bulk-referred Sh is taken over the step, [theta - dt, theta], with
Fo = 4 Dliq theta/dp^2, so the explicit source removes exactly the Newman
fraction of the step (newmanSphere::ShStep).
It is referred to the drop bulk, so use liquidDrive bulk; set donorSink drop
in the alpha.<gas> controls so the depletion layer stays sub-grid, and
dropHaloLayers (3) so smeared liquid around the drop ages and is depleted.
*/

#ifndef Foam_meltingEvaporationModels_gasDiffusionSimple_H
#define Foam_meltingEvaporationModels_gasDiffusionSimple_H

#include "InterfaceCompositionModel.H"
#include "volFields.H"

namespace Foam
{
class phasePair;

namespace meltingEvaporationModels
{

template<class Thermo, class OtherThermo>
class gasDiffusionSimple
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Interface reconstruction and dispersed-object geometry
        scalar isoAlpha_;
        scalar alphaBandMin_;
        scalar alphaGasSampleMin_;
        scalar alphaLiqSampleMin_;
        scalar alphaLo_;
        scalar alphaHi_;
        scalar alphaClusterCut_;
        label minDropCells_;
        label dropHaloLayers_;

        word dispersedPhaseName_;      
        scalar geometryVolumeFactor_;

    // Transfer-model constants
    dimensionedScalar C_;

    // Gas-side molecular diffusivity [m2/s]
    // Used for Sc_g, Sh_g and gas-side resistance.
    scalar Dgas_;

    // Liquid-side molecular diffusivity [m2/s]
    // Used for the liquid-side transfer resistance.
    scalar Dliq_;

    // Gas-side Sherwood closure
    word gasSherwoodModel_;
    scalar gasSherwood_;

    // Liquid-side Sherwood closure
    word liquidSherwoodModel_;
    scalar liquidSherwood_;

    // Liquid concentration in the driving force: interface (sampled) or bulk (drop mean)
    word liquidDrive_;

    // Gas viscosity, only needed for gas-side Ranz-Marshall
    scalar muGas_;

    scalar dp_;
    word dpMode_;

    // Per-drop transfer area: iso (iso-surface) or volumeSphere (pi*Dvol^2)
    word dropAreaMode_;

    // Rescale flux cells so each drop transfers through its full area
    bool correctDropArea_;

    // Drop slip velocity: interface (adjacent cells) or shell (drop - shell)
    word slipVelocityMode_;
    scalar slipShellInner_;
    scalar slipShellOuter_;

    // Contact age of the dispersed phase
    bool contactAge_;
    scalar ageAlphaCut_;
    word surfaceRenewal_;
    label ageTimeIndex_;
    autoPtr<volScalarField> contactAgePtr_;

    // Dimensionless Henry/partition ratio:
        scalar k_;

    // Optional diagnostics and output fields
        bool debugDeff_;
        bool debugDrops_;
        bool writeDropFields_;
        bool debugInternalCirculation_;


        volScalarField interfaceArea_;
        volScalarField mDotc_;
        volScalarField dropIdField_;
        volScalarField dropDpField_;
        volScalarField dropReField_;
        volScalarField dropShGasField_;
        volScalarField dropShLiqField_;

    // Private member functions
        void updateInterface(const volScalarField& T);

public:

    TypeName("gasDiffusionSimple");

    gasDiffusionSimple
    (
        const dictionary& dict,
        const phasePair& pair
    );

    virtual ~gasDiffusionSimple() = default;

    // Explicit mass transfer rate [kg/m3/s]
    virtual tmp<volScalarField> Kexp
    (
        const volScalarField& T
    ) override;

    // Zero linearized contribution; Kexp provides the explicit source.
    virtual tmp<volScalarField> KSp
    (
        label /*modelVariable*/,
        const volScalarField& field
    ) override
    {
        const fvMesh& mesh = this->mesh_;

        // Base: [kg/m3/s]
        dimensionSet dims = dimDensity/dimTime;

        // Detect what the caller is linearizing against by checking the passed field dims
        if (field.dimensions() == dimPressure)
        {
            dims /= dimPressure;
        }
        else if (field.dimensions() == dimTemperature)
        {
            dims /= dimTemperature;
        }

        return tmp<volScalarField>::New
        (
            IOobject
            (
                IOobject::groupName("KSp0", this->pair().name()),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dims, 0.0)
        );
    }

    // Zero explicit pressure/temperature linearization contribution.
    virtual tmp<volScalarField> KSu
    (
        label /*modelVariable*/,
        const volScalarField& /*field*/
    ) override
    {
        const fvMesh& mesh = this->mesh_;

        return tmp<volScalarField>::New
        (
            IOobject
            (
                IOobject::groupName("KSu0", this->pair().name()),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimDensity/dimTime, 0.0)
        );
    }

    // Latent heat (negligible here)
    virtual tmp<volScalarField> L
    (
        const word& /*speciesName*/,
        const volScalarField& /*T*/
    ) const override
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "L",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar(dimEnergy/dimMass, 0.0)
        );
    }

    virtual const dimensionedScalar& Tactivate() const noexcept override
    {
        static dimensionedScalar dummyT("T", dimTemperature, 0);
        return dummyT;
    }

    virtual bool includeDivU() const noexcept override
    {
        return false;
    }
};

} // namespace meltingEvaporationModels
} // namespace Foam

#ifdef NoRepository
#   include "gasDiffusionSimple.C"
#endif

#endif
// ************************************************************************* //
