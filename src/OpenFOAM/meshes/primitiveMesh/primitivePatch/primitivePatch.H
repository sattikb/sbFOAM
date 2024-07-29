//TYPEDEF FOR ADDRESSING A FACELIST SLICE

#ifndef primitivePatch_H
#define primitivePatch_H

#include "PrimitivePatch.H"
//#include "face.H"
//#include "SubList.H"
//#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef PrimitivePatch<SubList<face>, const pointField&> primitivePatch;
}

#endif
