/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "AhmedYovanovich.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AhmedYovanovich<CloudType>::AhmedYovanovich
(
    const dictionary& dict,
    CloudType& cloud
)
:
    HeatTransferModel<CloudType>(dict, cloud, typeName)
{}


template<class CloudType>
Foam::AhmedYovanovich<CloudType>::AhmedYovanovich(const AhmedYovanovich<CloudType>& htm)
:
    HeatTransferModel<CloudType>(htm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::AhmedYovanovich<CloudType>::Nu
(
    const scalar Re,
    const scalar Pr
) const
{
    const scalar gamma = min(1.0/pow(Re, 0.25), 1.0);

    return 2 + 0.775*pow(Re, 0.5)*pow(Pr, 1.0/3.0)/(sqrt(2*gamma + 1) * pow(1 + (1.0/(pow(2*gamma + 1, 3)*Pr)), 1.0/6.0));
}


// ************************************************************************* //
