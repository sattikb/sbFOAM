//  Distance calculation for cells with face on a wall.
//  Searches pointNeighbours to find closest.

#ifndef nearWallDistSB_H
#define nearWallDistSB_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class nearWallDist Declaration
\*---------------------------------------------------------------------------*/

class nearWallDist
:
    public volScalarField::Boundary
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Do all calculations
        void calculate();


public:

    // Constructors

        //- Construct from components
        nearWallDist(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        nearWallDist(const nearWallDist&) = delete;


    //- Destructor
    virtual ~nearWallDist();


    // Member Functions

        const volScalarField::Boundary& y() const
        {
            return *this;
        }

        //- Correct for mesh geom/topo changes
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nearWallDist&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
