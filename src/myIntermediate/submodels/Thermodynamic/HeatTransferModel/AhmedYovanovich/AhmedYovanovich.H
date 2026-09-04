/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::AhmedYovanovich

Group
    grpLagrangianIntermediateHeatTransferSubModels

Description
    Nusselt-number model using the empirical Ahmed-Yovanovich correlation
    to be used in modelling of the fluid-solid particle heat transfer coefficient:

    \f[
        \mathrm{Nu} = a + b \, \mathrm{Re}_p^{m} \, \mathrm{Pr}^{n}
        \mathrm{Nu} = 2 + \frac{0.775 \cdot \mathrm{Re}^{0.5} \cdot \mathrm{Pr}^{1/3}} {\sqrt{2\gamma + 1} \cdot \left[1 + \frac{1}{(2\gamma + 1)^3 \mathrm{Pr}}\right]^{1/6}}
    \f]
    with

    \f[
        \mathrm{Re}_p =
            \frac{\rho_c \, | \mathbf{u}_\mathrm{rel} | \, d_p}{\mu_c}
    \f]

    \f[
        \mathrm{Pr} = \frac{ C_p \, \mu_c }{ \kappa_c }
    \f]

    \f[
        \gamma = min(1, \frac{1}{\mathrm{Re}_p^{0.25}})
    \f]
    where

    \vartable
        \mathrm{Nu}    | Nusselt number
        \mathrm{Re}_p  | Particle Reynolds number
        \mathrm{Pr}    | Prandtl number
        d_p            | Particle diameter
        \rho_c         | Density of carrier in the film surrounding particle
        \mu_c          | Dynamic viscosity of carrier in the film surrounding particle
        \mathbf{u}_\mathrm{rel} | Relative velocity between particle and carrier
        \gamma         | Correlation coefficient
        C_p            | Specific heat capacity
        \kappa_c       | Thermal conductivity of carrier in the film
    \endvartable

    Reference:
    \verbatim
        Standard model:
            G. Refai Ahmed, M.M. Yovanovich, Approximate analytical solution of forced
            convection heat transfer from isothermal spheres for all Prandtl numbers, J.
            Heat Transfer 116 (4) (1994) 838–843
    \endverbatim

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        // Mandatory entries
        heatTransferModel    AhmedYovanovich;
    }
    \endverbatim

    where the entries mean:
    \table
      Property           | Description              | Type | Reqd | Deflt
      heatTransferModel  | Type name: AhmedYovanovich  | word | yes  | -
    \endtable

SourceFiles
    AhmedYovanovich.C

\*---------------------------------------------------------------------------*/

#ifndef AhmedYovanovich_H
#define AhmedYovanovich_H

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class AhmedYovanovich Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class AhmedYovanovich
:
    public HeatTransferModel<CloudType>
{
    // Private Data


public:

    //- Runtime type information
    TypeName("AhmedYovanovich");


    // Generated Methods

        //- No copy assignment
        void operator=(const AhmedYovanovich&) = delete;


    // Constructors

        //- Construct from dictionary
        AhmedYovanovich(const dictionary& dict, CloudType& cloud);

        //- Copy construct
        AhmedYovanovich(const AhmedYovanovich<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<HeatTransferModel<CloudType>> clone() const
        {
            return autoPtr<HeatTransferModel<CloudType>>
            (
                new AhmedYovanovich<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~AhmedYovanovich() = default;


    // Member Functions

        // Evaluation

            //- Return Nusselt number
            virtual scalar Nu
            (
                const scalar Re,
                const scalar Pr
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AhmedYovanovich.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
