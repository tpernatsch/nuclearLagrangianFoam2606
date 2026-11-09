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
    liquidSherwoodModel     Higbie;
    liquidSherwood          2;

    dpMode                  perDropVA;
    dp                      3.0e-4;
}
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

    // Gas viscosity, only needed for gas-side Ranz-Marshall
    scalar muGas_;

    scalar dp_;
    word dpMode_;

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
