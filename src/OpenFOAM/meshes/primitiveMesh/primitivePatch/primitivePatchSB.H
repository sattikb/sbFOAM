//TYPEDEF FOR ADDRESSING A FACELIST SLICE

#ifndef primitivePatchSB_H
#define primitivePatchSB_H

#include "PrimitivePatch.H"

namespace Foam
{
    typedef PrimitivePatch<SubList<face>, const pointField&> primitivePatch;
}

#endif
