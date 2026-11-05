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
    Foam::EkanayakeLiftForce

Group
    grpLagrangianIntermediateForceSubModels

Description
    Shear and wall-induced lift on a small sphere near a wall, valid for any
    particle-wall distance (inner and outer region), from Ekanayake, Berry &
    Harvie (2021), Appendix A, freely rotating sphere.

    With a = d/2, the slip along the wall u_s = (U_p - U_c) - normal part,
    the wall shear gamma = dU_t/dn (n into the fluid, t along u_s),
    Re_slip = |u_s| a/nu, Re_gamma = |gamma| a^2/nu, l* = l/a
    (l = distance of the parcel centre from the nearest wall face),
    l/L_G = l* sqrt(Re_gamma), l/L_S = l* Re_slip,
    epsilon = sqrt(Re_gamma)/Re_slip:
    \f[
        F_L \frac{\rho}{\mu^2} = C_{L,1} Re_\gamma^2
          + sgn(\gamma) C_{L,2} Re_{slip} Re_\gamma
          + C_{L,3} Re_{slip}^2        \qquad (+ = away from the wall)
    \f]
    The outer part of C_{L,2}, -f(eps, l/L_G) (9/pi) eps J(eps)/|gamma*|,
    is f times the unbounded McLaughlin lift; it is applied here in its
    vector form (as McLaughlinMeiLiftForce, F ~ (U_c - U_p) x curl U_c), so
    that away from walls (f -> 1, C_{L,1}, C_{L,3} -> 0) the force reduces
    to the unbounded shear lift in any flow. The rest acts along n:
    \f[
        C_{L,1} = f_1 \, 1.982\,e^{-0.115 x^2 - 0.2771 x}
                + f_2 (0.89934/l^* - 1.9610/l^{*2} + 1.0216/l^{*3})
    \f]
    with x = l/L_G, f_1 = 0.925 e^{-0.35 Re_gamma} - 0.0135 e^{-7000
    Re_gamma}, f_2 = 1 + sqrt(Re_gamma) (Ekanayake et al. 2020, eqs 4.1-4.5,
    Table 3, freely rotating; note: Appendix A of the 2021 paper lists the
    non-rotating lambda_8, lambda_9 by mistake),
    \f[
        C_{L,2,in} = -2.6729 - 0.8373/l^* + 0.4683/l^{*2}
    \f]
    (Cherukat & McLaughlin 1995, without the l* term, which the outer part
    replaces),
    \f[
        C_{L,3} = \frac{18\pi}{32 + 2 y + 3.8 y^2 + 0.049 y^3}
                + 0.4757/l^* - 1.268/l^{*2} + 0.683/l^{*3}
    \f]
    with y = l/L_S (Takemura 2004; Ekanayake et al. 2021 eq 4.4b), and
    \f[
        J(\epsilon) = 2.255 (1 + 0.2\epsilon^{-2})^{-3/2},
        \qquad
        f(\epsilon, x) = 1 - \exp[-\tfrac{11}{96}\pi^2 x/J(\epsilon)]
    \f]
    (Legendre & Magnaudet 1998; Takemura, Magnaudet & Dimitrakopoulos 2009).

    hinderedMigration: the drag on a sphere moving towards a wall is larger
    than the (unbounded) drag of the drag model. Since the parcels move at
    their terminal velocity (dt >> tau_p), the wall-normal part of the lift
    is scaled by 6 pi/C_D,perp (Brenner; Ekanayake et al. 2021 eq 5.4):
    \f[
        h_\perp = 1 - \tfrac{9}{8} s + \tfrac{1}{2} s^3
                - \tfrac{135}{256} s^4 - \tfrac{1}{8} s^5,
        \qquad s = 1/l^*
    \f]
    giving the hindered drift w = F_L/(mu a C_D,perp) of the paper.

    slip: "relative" (default) uses U_p - U_c. With dt >> tau_p, U_p is the
    end-of-step velocity of the previous step while U_c is taken at the new
    position, so for small heavy particles the change of U_c over one step
    is of the order of the slip itself (10 um in the heated pool: 10-90 %
    spread as large as the mean). "terminal" uses the settling velocity
    relative to the fluid, u_t g/|g| with (rho_p - rho_c) g d^2/(18 mu_c)
    = u_t (1 + 0.15 Re^0.687) (Schiller-Naumann, as sphereDrag), i.e. the
    slip of a parcel at its terminal state; it neglects the inertial lag
    tau_p DU_c/Dt (a few % for 100 um here) and slip from other forces.
    Ekanayake et al. (2021, eq 5.2) obtain the slip the same way when
    comparing with experiments.

    Finite Reynolds number: the correlation above is fitted for Re_slip,
    Re_gamma <= 0.1. When either exceeds ReMax (default 0.1), the
    decomposition of Takemura, Magnaudet & Dimitrakopoulos (2009) with the
    correlations recommended by Shi & Rzehak (2020, eq 24, Table 3) is used
    instead (diameter-based Re_p = 2 Re_slip, L_u = l/L_S, L_omega = l/L_G,
    L_R = l*):
    \f[
        F_L = C^W_{L\omega} f_L \tfrac{\pi}{8} d^2 \rho_c |u_s|^2 \, n
            + h_L F_{ub}
    \f]
    with the wall lift of Zeng et al. (2009) on the low-Re base of Shi &
    Rzehak (eq 20),
    \f[
        C^W_{L\omega} = c(L_R) + [C_{1.01} - c(1.01)]
            e^{-11 (\delta/\delta_c)^{1.2}},
        \quad c(L_R) = C^W_0(L_u) A^{(L_R/2)^g}
    \f]
    C^W_0 = 9/8/(1 + 0.13 L_u (L_u + 0.53)) (L_u < 10), 8.94 L_u^-2.09
    (L_u >= 10), A = 1 + 0.329 Re_p + 0.00485 Re_p^2,
    g = -0.9 tanh(0.022 Re_p), delta = (L_R - 1)/2,
    delta_c = 3 exp(-0.17 Re_p^0.7), C_1.01 = 0.313 + 0.812 exp(-0.125
    Re_p^0.77); the shear modification of the wall lift
    f_L = exp(-0.22 eps^0.8 L_omega^2.5) (Takemura et al. 2009, eq 30 of
    Shi & Rzehak) and the wall modification of the shear lift
    h_L = 1 - exp[-11/96 pi^2 70/(70 + L_omega^1.378) L_omega/J]
    (L_omega < 15), 1 - 1.8778 L_omega^(-5/3)/J (L_omega >= 15)
    (Shi & Rzehak eq 29), times the unbounded lift F_ub (unbounded shear
    lift 18/pi^2 eps J, Legendre & Magnaudet). Shi & Rzehak find it
    reasonable up to Re_p ~ 4 (Re_slip ~ 2); a warning is printed once above.

    Limits: l* >= 1.2 (l* is clamped there); a single plane wall (the
    nearest one is used); slip parallel to the wall. Parcels smaller than
    dMin get no lift (also keeps Brownian kicks out of the slip).

    scale: multiplies the whole force, shear and wall parts (default 1). It
    has no physical meaning: it is for testing the sensitivity of a result
    to the uncertainty of the lift.

    References:
        Ekanayake, N.I.K., Berry, J.D., Harvie, D.J.E. (2021). Lift and drag
        forces acting on a particle moving in the presence of slip and shear
        near a wall. J. Fluid Mech. 915, A103.
        Ekanayake, N.I.K., Berry, J.D., Stickland, A.D., Dunstan, D.E.,
        Muir, I.L., Dower, S.K., Harvie, D.J.E. (2020). Lift and drag forces
        acting on a particle moving with zero slip in a linear shear flow
        near a wall. J. Fluid Mech. 904, A6.
        Shi, P., Rzehak, R. (2020). Lift forces on solid spherical particles
        in wall-bounded flows. Chem. Eng. Sci. 211, 115264.
        Takemura, F., Magnaudet, J., Dimitrakopoulos, P. (2009). Migration
        and deformation of bubbles rising in a wall-bounded shear flow at
        finite Reynolds number. J. Fluid Mech. 634, 463-486.
        Zeng, L., Najjar, F., Balachandar, S., Fischer, P. (2009). Forces on
        a finite-sized particle located close to a wall in a linear shear
        flow. Phys. Fluids 21, 033302.

Usage
    \verbatim
    particleForces
    {
        EkanayakeLiftForce
        {
            dMin                1e-6;       // [m] no lift below (default 0)
            hinderedMigration   yes;        // default yes
            slip                terminal;   // relative (default) | terminal
            ReMax               0.1;        // above: Takemura et al. (2009)
            scale               1;          // default 1, sensitivity tests
            patches             (innerWall outerWall bottom); // default: walls
        }
    }

    interpolationSchemes
    {
        curlUcDt        cellPatchConstrained;
        gradUcDt        cellPatchConstrained;
    }
    \endverbatim

SourceFiles
    EkanayakeLiftForce.C

\*---------------------------------------------------------------------------*/

#ifndef EkanayakeLiftForce_H
#define EkanayakeLiftForce_H

#include "LiftForce.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class EkanayakeLiftForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EkanayakeLiftForce
:
    public LiftForce<CloudType>
{
    // Private Data

        //- No lift below this diameter [m]
        const scalar dMin_;

        //- Scale the wall-normal lift by the wall-normal drag hindrance
        const bool hinderedMigration_;

        //- Re_slip or Re_gamma above which the Takemura et al. (2009)
        //  decomposition replaces the Ekanayake correlation
        const scalar ReMax_;

        //- Multiplier of the whole force (sensitivity tests, default 1)
        const scalar scale_;

        //- Wall patches (empty: all wall patches)
        const wordRes patches_;

        //- Slip = terminal settling velocity instead of U_p - U_c
        bool terminalSlip_;

        //- Unit normal (into the fluid) of the nearest wall face, per cell;
        //  zero where no wall is reached
        vectorField nWall_;

        //- Centre of the nearest wall face, per cell
        vectorField CfWall_;

        //- Gradient of the carrier velocity interpolator
        autoPtr<interpolation<tensor>> gradUcInterpPtr_;

        //- Out-of-range warning (Re_p > 4) printed
        mutable bool warned_;


    // Private Member Functions

        //- Nearest-wall normal and face centre for every cell
        void setWallData();

        //- Terminal settling velocity relative to the fluid (along g)
        vector terminalSlip
        (
            const scalar d,
            const scalar rhop,
            const scalar rhoc,
            const scalar muc
        ) const;

        //- Legendre & Magnaudet (1998) J(epsilon)
        static scalar J(const scalar epsilon);

        //- Wall lift coefficient at Re_p << 1 vs L_u (Shi & Rzehak eq 20)
        static scalar CW0(const scalar Lu);

        //- Wall lift coefficient at finite Re_p (Zeng et al. 2009)
        static scalar CW(const scalar LR, const scalar Rep);


public:

    //- Runtime type information
    TypeName("EkanayakeLiftForce");


    // Constructors

        //- Construct from mesh
        EkanayakeLiftForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& forceType = typeName
        );

        //- Construct copy
        EkanayakeLiftForce(const EkanayakeLiftForce& lf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new EkanayakeLiftForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~EkanayakeLiftForce() = default;


    // Member Functions

        //- Cache fields
        virtual void cacheFields(const bool store);

        //- Calculate the coupled force
        virtual forceSuSp calcCoupled
        (
            const typename CloudType::parcelType& p,
            const typename CloudType::parcelType::trackingData& td,
            const scalar dt,
            const scalar mass,
            const scalar Re,
            const scalar muc
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EkanayakeLiftForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
