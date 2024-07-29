// CONTAINER CLASSES FOR TYPE BOOLEAN LISTS

#ifndef boolList_H
#define boolList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef UList<bool> boolUList;

    typedef List<bool> boolList;
    typedef List<List<bool>> boolListList;
}

#endif
