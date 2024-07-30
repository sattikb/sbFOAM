// Line using referred points.

#ifndef linePointRefSB_H
#define linePointRefSB_H

#include "lineSB.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef line<point, const point&> linePointRef;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
