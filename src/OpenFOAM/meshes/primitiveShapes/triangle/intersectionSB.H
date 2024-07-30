// PRIMARILY TO DETERMINE IF A HALF-RAY FROM A POINT IN PLANE OF TRIANGLE
// INTERSECTS THE TRIANGLE OR NOT.

#ifndef intersectionSB_H
#define intersectionSB_H

#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class intersection Declaration
\*---------------------------------------------------------------------------*/

class intersection
{
    // Static data

        //- Relative tolerance for point in a plane.
        //  The tolerance is relative to the object size.
        //  Used to determine if a half-ray from a point in plane of triangle
        // intersects the triangle or not.
        static scalar planarTol_;

public:

    enum class direction
    {
        vector,
        contactSphere
    };

    enum class algorithm
    {
        fullRay,        // Intersecting with full ray
        halfRay,        // Intersecting with half ray
        visible         // Intersecting with the visible side
    };

    // Static Member Functions

        //- Direction names
        static const NamedEnum<direction, 2> directionNames_;

        //- Projection algorithm names
        static const NamedEnum<algorithm, 3> algorithmNames_;

        //- Return planar tolerance
        static scalar planarTol()
        {
            return planarTol_;
        }

        //- Set the planar tolerance, returning the previous value
        static scalar setPlanarTol(const scalar t)
        {
            if (t < -vSmall)
            {
                FatalErrorInFunction << abort(FatalError);
            }

            scalar oldTol = planarTol_;
            planarTol_ = t;

            return oldTol;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
