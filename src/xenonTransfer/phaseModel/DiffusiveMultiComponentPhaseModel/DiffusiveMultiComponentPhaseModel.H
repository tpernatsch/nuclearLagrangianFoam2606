/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::DiffusiveMultiComponentPhaseModel

Description
    Phase model with multiple species. Provides access to species mass
    fractions and solves their governing equations.

SourceFiles
    DiffusiveMultiComponentPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef DiffusiveMultiComponentPhaseModel_H
#define DiffusiveMultiComponentPhaseModel_H

#include "phaseModel.H"
#include "hashedWordList.H"
#include "rhoReactionThermo.H"
#include "basicThermo.H"

namespace Foam
{
class multiphaseInterSystem;

template<class BasePhaseModel, class phaseThermo>
class DiffusiveMultiComponentPhaseModel
:
    public BasePhaseModel
{
protected:

    //- Species table
    hashedWordList species_;

    //- Inert species index
    label inertIndex_;

    //- Thermophysical model
    autoPtr<phaseThermo> thermoPtr_;

    //- Volume fractions for species
    PtrList<volScalarField> X_;

    //- Add diffusion term to Yi equations
    bool addDiffusion_;

    //- Molecular diffusivity (m2/s)
    dimensionedScalar Dmol_;

    //- Turbulent Schmidt number
    scalar Sct_;

    //- Conservative species storage, Z = alpha*Y
    PtrList<volScalarField> Z_;

    //- Transform volume fractions into mass fractions
    void calculateMassFractions();

    //- Transform mass fractions into volume fractions
    void calculateVolumeFractions();

    //- Synchronise one Z boundary from alpha*Y
    void syncZBoundary(const label i);

    //- Synchronise one full Z field from alpha*Y
    void syncZField(const label i);


public:

    // Constructors
    DiffusiveMultiComponentPhaseModel
    (
        const multiphaseInterSystem& fluid,
        const word& phaseName
    );

    //- Destructor
    virtual ~DiffusiveMultiComponentPhaseModel() = default;


    // Access

    const hashedWordList& species() const
    {
        return species_;
    }


    // Thermo

    virtual const phaseThermo& thermo() const;
    virtual phaseThermo& thermo();

    virtual void correct();

    virtual void solveYi
    (
        PtrList<volScalarField::Internal>& Su,
        PtrList<volScalarField::Internal>& Sp
    );

    //- Make Z/Y consistent with the current alpha before Kexp sampling
    bool preKexpSyncSpeciesStorage
    (
        const word& speciesName,
        const bool verbose = false
    );

    virtual const PtrList<volScalarField>& Y() const;
    virtual PtrList<volScalarField>& Y();

    label inertIndex() const;
};

} // End namespace Foam


#ifdef NoRepository
#include "DiffusiveMultiComponentPhaseModel.C"
#endif

#endif

// ************************************************************************* //
