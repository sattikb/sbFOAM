/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Primitive
    bool

Description
    System bool

SourceFiles
    boolIO.C

\*---------------------------------------------------------------------------*/

#ifndef bool_H
#define bool_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Istream& operator>>(Istream&, bool&);
Ostream& operator<<(Ostream&, const bool);

bool readBool(Istream&);

inline void writeEntry(Ostream& os, const bool value)
{
    os << value;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pTraits.H"
#include "directionSB.H"

namespace Foam
{

// template specialisation for pTraits<bool>
template<>
class pTraits<bool>
{
    bool p_;

public:

    //- Component type
    typedef bool cmptType;

    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of bool is 0
        static const direction rank = 0;

        //- Number of components in bool is 1
        static const direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const bool zero;
        static const bool one;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const bool&);

        //- Construct from Istream
        pTraits(Istream&);


    // Member Functions

        //- Access to the bool value
        operator bool() const
        {
            return p_;
        }

        //- Access to the bool value
        operator bool&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
