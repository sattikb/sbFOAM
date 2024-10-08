/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Namespace
    Foam::SRF

Description
    Namespace for single rotating frame (SRF) models

Class
    Foam::SRF::SRFModel

Description
    Top level model for single rotating frame
    - Steady state only - no time derivatives included

SourceFiles
    SRFModel.C

\*---------------------------------------------------------------------------*/

#ifndef SRFModel_H
#define SRFModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volFields.H"
#include "vectorFieldSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SRF
{

/*---------------------------------------------------------------------------*\
                          Class SRFModel Declaration
\*---------------------------------------------------------------------------*/

class SRFModel
:
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to the relative velocity field
        const volVectorField& Urel_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Origin of the axis
        dimensionedVector origin_;

        //- Axis of rotation, a direction vector which passes through the origin
        vector axis_;

        //- SRF model coefficients dictionary
        dictionary SRFModelCoeffs_;

        //- Angular velocity of the frame (rad/s)
        dimensionedVector omega_;


public:

    //- Runtime type information
    TypeName("SRFModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             SRFModel,
             dictionary,
             (
                 const volVectorField& Urel
             ),
             (Urel)
         );


    // Constructors

        //- Construct from components
        SRFModel
        (
            const word& type,
            const volVectorField& Urel
        );

        //- Disallow default bitwise copy construction
        SRFModel(const SRFModel&) = delete;


    // Selectors

         //- Return a reference to the selected SRF model
         static autoPtr<SRFModel> New
         (
             const volVectorField& Urel
         );


    //- Destructor
    virtual ~SRFModel();


    // Member Functions

        // Edit

            //- Read SRFProperties dictionary
            virtual bool read();


        // Access

            //- Return the origin of rotation
            const dimensionedVector& origin() const;

            //- Return the axis of rotation
            const vector& axis() const;

            //- Return the angular velocity field [rad/s]
            const dimensionedVector& omega() const;

            //- Return the coriolis force
            tmp<volVectorField::Internal> Fcoriolis() const;

            //- Return the centrifugal force
            tmp<volVectorField::Internal> Fcentrifugal() const;

            //- Source term component for momentum equation
            tmp<volVectorField::Internal> Su() const;

            //- Return velocity vector from positions
            vectorField velocity(const vectorField& positions) const;

            //- Return velocity of SRF for complete mesh
            tmp<volVectorField> U() const;

            //- Return absolute velocity for complete mesh
            tmp<volVectorField> Uabs() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SRFModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SRF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
