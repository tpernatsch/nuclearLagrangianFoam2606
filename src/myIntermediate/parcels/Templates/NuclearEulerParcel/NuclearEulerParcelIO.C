/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tommaso Pernatsch
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

Description
    NuclearEulerParcel adds exactly one new persisted per-parcel field
    relative to its ParcelType (NuclearParcel): the phase-switch parameter
    chi_. Every other extra quantity it reads (the secondary Euler phase's
    fields, and the primary phase's turbulence fields) stays transient,
    carried only in trackingData and recomputed by interpolation every
    tracking sub-step - only chi_ needs real I-O.

\*---------------------------------------------------------------------------*/

#include "NuclearEulerParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::NuclearEulerParcel<ParcelType>::propertyList_ =
    Foam::NuclearEulerParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::NuclearEulerParcel<ParcelType>::sizeofFields
(
    sizeof(NuclearEulerParcel<ParcelType>)
  - offsetof(NuclearEulerParcel<ParcelType>, chi_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::NuclearEulerParcel<ParcelType>::NuclearEulerParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    chi_(false)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is >> chi_;
        }
        else
        {
            is.fatalCheckNativeSizes(FUNCTION_NAME);

            is.read(reinterpret_cast<char*>(&chi_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::NuclearEulerParcel<ParcelType>::readFields(CloudType& c)
{
    const bool readOnProc = c.size();

    ParcelType::readFields(c);

    IOField<bool> chi(c.fieldIOobject("chi", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, chi);

    label i = 0;
    for (NuclearEulerParcel<ParcelType>& p : c)
    {
        p.chi_ = chi[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::NuclearEulerParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<bool> chi(c.fieldIOobject("chi", IOobject::NO_READ), np);

    label i = 0;
    for (const NuclearEulerParcel<ParcelType>& p : c)
    {
        chi[i] = p.chi_;

        ++i;
    }

    chi.write(writeOnProc);
}


template<class ParcelType>
void Foam::NuclearEulerParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                          \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("chi", chi_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::NuclearEulerParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, obr);

    if (!c.size()) return;

    auto& chi = cloud::lookupIOField<bool>("chi", obr);

    label i = 0;
    for (NuclearEulerParcel<ParcelType>& p : c)
    {
        p.chi_ = chi[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::NuclearEulerParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    auto& chi = cloud::createIOField<bool>("chi", np, obr);

    label i = 0;
    for (const NuclearEulerParcel<ParcelType>& p : c)
    {
        chi[i] = p.chi_;

        ++i;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const NuclearEulerParcel<ParcelType>& p
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.chi();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.chi_),
            NuclearEulerParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
