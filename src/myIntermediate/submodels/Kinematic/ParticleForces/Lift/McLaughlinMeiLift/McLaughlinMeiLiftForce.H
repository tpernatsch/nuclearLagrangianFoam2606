/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::McLaughlinMeiLiftForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Shear lift on a small sphere with McLaughlin's correction for a slip
    that is not small compared with the shear (Saffman's theory needs
    epsilon >> 1), using Mei's fit of McLaughlin's results.

    Saffman (1965, 1968):
    \f[
        F_S = 6.46\,\mu\,a^2\,|U_s|\,\sqrt{|\omega|/\nu},
        \qquad a = d/2
    \f]
    McLaughlin (1991):
    \f[
        F_L = F_S\,\frac{J(\epsilon)}{J_\infty},
        \qquad
        \epsilon = \frac{\sqrt{|\omega|\nu}}{|U_s|}
                 = \frac{\sqrt{Re_\omega}}{Re},
        \qquad J_\infty = 2.255
    \f]
    Mei (1992), fit of McLaughlin's J for 0.1 <= epsilon <= 20:
    \f[
        J(\epsilon) = 0.3\,J_\infty
            \left\{1 + \tanh\left[\tfrac{5}{2}
                \left(\log_{10}\epsilon + 0.191\right)\right]\right\}
            \left\{\tfrac{2}{3} + \tanh\left(6\epsilon - 1.92\right)\right\}
    \f]
    with Re = rho_c |U_c - U_p| d / mu_c (slip Reynolds number) and
    Re_omega = rho_c |curl U_c| d^2 / mu_c (shear Reynolds number).
    J -> J_inf (Saffman) for large epsilon; J < 0 (lift reversed) for
    epsilon below about 0.19, where the fit tends to zero.

    In the form used by LiftForce, F = (m/rho_p) rho_c Cl (U_c - U_p) x
    curl(U_c), this gives
    \f[
        Cl = \frac{27}{2\pi^2}\,\frac{J(\epsilon)}{\sqrt{Re_\omega}}
    \f]
    (= the SaffmanMeiLiftForce coefficient 3 x 6.46/(2 pi sqrt(Re_omega))
    at J = J_inf, without Mei's finite-Re factor).

    McLaughlin's theory is for small Re and Re_omega (here < 1); for Re of
    order 10 and above use SaffmanMeiLiftForce.

    References:
        Saffman, P.G. (1965). The lift on a small sphere in a slow shear
        flow. J. Fluid Mech. 22, 385-400; corrigendum (1968) 31, 624.
        McLaughlin, J.B. (1991). Inertial migration of a small sphere in
        linear shear flows. J. Fluid Mech. 224, 261-274.
        Mei, R. (1992). An approximate expression for the shear lift force
        on a spherical particle at finite Reynolds number. Int. J.
        Multiphase Flow 18(1), 145-147.

Usage
    \verbatim
    particleForces
    {
        McLaughlinMeiLiftForce {}
    }
    \endverbatim

SourceFiles
    McLaughlinMeiLiftForce.C

\*---------------------------------------------------------------------------*/

#ifndef McLaughlinMeiLiftForce_H
#define McLaughlinMeiLiftForce_H

#include "LiftForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class McLaughlinMeiLiftForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class McLaughlinMeiLiftForce
:
    public LiftForce<CloudType>
{
protected:

    // Protected Member Functions

        //- McLaughlin's J(epsilon), Mei's fit
        static scalar J(const scalar epsilon);

        //- Calculate the lift coefficient
        virtual scalar Cl
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const vector& curlUc,
            const scalar Re,
            const scalar muc
        ) const;


public:

    //- Runtime type information
    TypeName("McLaughlinMeiLiftForce");


    // Constructors

        //- Construct from mesh
        McLaughlinMeiLiftForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        McLaughlinMeiLiftForce(const McLaughlinMeiLiftForce& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new McLaughlinMeiLiftForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~McLaughlinMeiLiftForce() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "McLaughlinMeiLiftForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
