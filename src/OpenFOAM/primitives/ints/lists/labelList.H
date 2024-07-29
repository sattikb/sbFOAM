// PREVIOUSLY TEMPLATED LIST, NOW DEFINED FOR LABELS

#ifndef labelList_H
#define labelList_H

//#include "label.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Note: frequently used UList version is located in container itself

    typedef List<label> labelList;
    typedef List<labelList> labelListList;
    typedef List<labelListList> labelListListList;
}

#endif
