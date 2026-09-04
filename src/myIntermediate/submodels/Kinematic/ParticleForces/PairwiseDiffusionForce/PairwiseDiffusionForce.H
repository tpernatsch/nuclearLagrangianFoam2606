/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::PairWiseDiffusionForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Particle diffusion force model involving pairwise interactions, using theory from smoothed particle hydrodynamics (SPH).

    \f[
        \vec{F}_\mathrm{Dif} = \sum_neighbors D \cdot (\vec{x_j} - \vec{x_i}) \cdot W(r_{ij}, h)
    \f]

    where

    \vartable
        \vec{f}_\mathrm{Dif}  | Specific diffusion force       [m s^{-2}]
        D                     | Diffusivity coefficient        [s^{-2}]
        \vec{x_i}             | Position of particle i         [m]
        \vec{x_j}             | Position of particle j         [m]
        W(r_{ij}, h)          | Kernel function             [-]
    \endvartable

Usage
    Minimal example by using \c constant/cloudProperties:
    \verbatim
    subModels
    {
        particleForces
        {
            pairwiseDiffusion
            {
                type    pairwiseDiffusion;
                D       <scalar>;
                h       <scalar>;
            }
        }
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                       | Type   | Reqd | Deflt
      type     | Type name: pairwiseDiffusion      | word   | yes  | -
      D        | Diffusivity coefficient           | scalar | yes  | -
      h        | Cutoff radius                     | scalar | yes  | -
    \endtable

SourceFiles
    PairwiseDiffusionForce.C

\*---------------------------------------------------------------------------*/

#ifndef PairwiseDiffusionForce_H
#define PairwiseDiffusionForce_H

#include "ParticleForce.H"

namespace Foam
{

template<class CloudType>
class PairwiseDiffusionForce : public ParticleForce<CloudType>
{
    // Private data
    const scalar D_;              // Diffusivity
    const scalar h_;              // Cutoff radius

public:

    TypeName("pairwiseDiffusion");

    // Constructors
    PairwiseDiffusionForce(CloudType& owner, const fvMesh& mesh, const dictionary& dict);

    PairwiseDiffusionForce(const PairwiseDiffusionForce& dff);

    virtual forceSuSp calcNonCoupled
    (
        const typename CloudType::parcelType& p,
        const typename CloudType::parcelType::trackingData& td,
        const scalar dt,
        const scalar mass,
        const scalar Re,
        const scalar muc
    ) const;
};

} // namespace Foam

#ifdef NoRepository
    #include "PairwiseDiffusionForce.C"
#endif

#endif


// ************************************************************************* //
