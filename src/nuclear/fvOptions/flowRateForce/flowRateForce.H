/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fv::flowRateForce

Description
    Calculates and applies the force necessary to maintain the specified
    flow rate across the specified faceZone.

    Note: Currently only handles kinematic pressure (incompressible solvers).

Usage
    Example usage:
    \verbatim
    selectionMode   all;                    // Apply force to all cells
    fields          (U);                    // Name of velocity field
    direction       (1 0 0);                // Direction of the applied source
    faceZone        facesZoneInlet          // Name of the faceZone
    flowRate        1;                      // Value of the flow rate
    relaxation      0.1;                    // Optional relaxation factor
    \endverbatim


Note
- "direction" is used to specify the direction vector of the momentum source;
  the provided value is normalized to ensure unit magnitude
- "flowRate" must be specified in m^3/s
- the absolute value of the calculated faceZone flow rate is used;
  no particular care about face orientation should be needed


SourceFiles
    flowRateForce.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateForce_H
#define flowRateForce_H

#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class flowRateForce Declaration
\*---------------------------------------------------------------------------*/

class flowRateForce
:
    public cellSetOption
{
protected:

    // Protected data

        //- Flow rate
        scalar G_;

        //- Pressure gradient before correction
        scalar gradP0_;

        //- Change in pressure gradient
        scalar dGradP_;

        //- Flow vector
        vector flowVec_;

        //- Flow direction
        vector flowDir_;

        //- Relaxation factor
        scalar relaxation_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> rAPtr_;

        //- Name of the faceZone
        word faceZoneName_;

        //- Id for the face zone
        label zoneID_;

         //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- List of +1/-1 representing face flip map (1 use as is, -1 negate)
        labelList faceSign_;

        //- Area of the face zone
        scalar faceZoneArea_;


    // Protected Member Functions

        //- Calculate and return the absolute value of the flow rate
        //  integrated over the selected faceZone
        virtual scalar magFlowRate() const;

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        flowRateForce(const flowRateForce&);

        //- Disallow default bitwise assignment
        void operator=(const flowRateForce&);

        //- Initialise source model
        void initialise();

        //- Calculate total area of faceZone across processors
        void calculateTotalArea(scalar& area);


public:

    //- Runtime type information
    TypeName("flowRateForce");


    // Constructors

        //- Construct from explicit source name and mesh
        flowRateForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
