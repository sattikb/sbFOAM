/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
    Foam::AMIMethod

Description
    Base class for Arbitrary Mesh Interface (AMI) methods

SourceFiles
    AMIMethod.C

\*---------------------------------------------------------------------------*/

#ifndef AMIMethod_H
#define AMIMethod_H

#include "className.H"
#include "DynamicList.H"
#include "faceAreaIntersect.H"
#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"
#include "treeBoundBoxList.H"
#include "primitivePatchSB.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class AMIMethod Declaration
\*---------------------------------------------------------------------------*/

class AMIMethod
{
protected:

    //- Local typedef to octree tree-type
    typedef treeDataPrimitivePatch<primitivePatch> treeType;


    // Protected data

        //- Reference to source patch
        const primitivePatch& srcPatch_;

        //- Reference to target patch
        const primitivePatch& tgtPatch_;

        //- Flag to indicate that the two patches are co-directional and
        //  that the orientation of the target patch should be reversed
        const bool reverseTarget_;

        //- Flag to indicate that the two patches must be matched/an overlap
        //  exists between them
        const bool requireMatch_;

        //- Source face areas
        const scalarField& srcMagSf_;

        //- Target face areas
        const scalarField& tgtMagSf_;

        //- Labels of faces that are not overlapped by any target faces
        //  (should be empty for correct functioning)
        labelList srcNonOverlap_;

        //- Octree used to find face seeds
        autoPtr<indexedOctree<treeType>> treePtr_;

        //- Face triangulation mode
        const faceAreaIntersect::triangulationMode triMode_;


    // Protected Member Functions

        // Helper functions

            //- Check AMI patch coupling
            void checkPatches() const;

            //- Initialise and return true if all ok
            bool initialise
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label& srcFacei,
                label& tgtFacei
            );

            //- Write triangle intersection to OBJ file
            void writeIntersectionOBJ
            (
                const scalar area,
                const face& f1,
                const face& f2,
                const pointField& f1Points,
                const pointField& f2Points
            ) const;


        // Common AMI method functions

            //- Reset the octree for the target patch face search
            void resetTree();

            //- Find face on target patch that overlaps source face
            label findTargetFace(const label srcFacei) const;

            //- Add faces neighbouring facei to the ID list
            void appendNbrFaces
            (
                const label facei,
                const primitivePatch& patch,
                const DynamicList<label>& visitedFaces,
                DynamicList<label>& faceIDs
            ) const;

            //- The maximum edge angle that the walk will cross
            virtual scalar maxWalkAngle() const;


public:

    //- Runtime type information
    TypeName("AMIMethod");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        AMIMethod,
        components,
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget,
            const bool requireMatch
        ),
        (
            srcPatch,
            tgtPatch,
            srcMagSf,
            tgtMagSf,
            triMode,
            reverseTarget,
            requireMatch
        )
    );


    // Constructors

        //- Construct from components
        AMIMethod
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const scalarField& srcMagSf,
            const scalarField& tgtMagSf,
            const faceAreaIntersect::triangulationMode& triMode,
            const bool reverseTarget,
            const bool requireMatch
        );

        //- Disallow default bitwise copy construction
        AMIMethod(const AMIMethod&) = delete;


    //- Selector
    static autoPtr<AMIMethod> New
    (
        const word& methodName,
        const primitivePatch& srcPatch,
        const primitivePatch& tgtPatch,
        const scalarField& srcMagSf,
        const scalarField& tgtMagSf,
        const faceAreaIntersect::triangulationMode& triMode,
        const bool reverseTarget,
        const bool requireMatch
    );


    //- Destructor
    virtual ~AMIMethod();


    // Member Functions

        // Access

            //- Labels of faces that are not overlapped by any target faces
            //  Note: this should be empty for correct functioning
            inline const labelList& srcNonOverlap() const;

            //- Flag to indicate that interpolation patches are conformal
            virtual bool conformal() const;


        // Manipulation

            //- Update addressing and weights
            virtual void calculate
            (
                labelListList& srcAddress,
                scalarListList& srcWeights,
                labelListList& tgtAddress,
                scalarListList& tgtWeights,
                label srcFacei = -1,
                label tgtFacei = -1
            ) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const AMIMethod&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "AMIMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
