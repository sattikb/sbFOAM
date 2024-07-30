// TYPEF OF THE 1-D LISTS CONTAINING INFORMATION ON THE EDGES

#ifndef edgeList_H
#define edgeList_H

#include "edgeSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef List<edge> edgeList;
    typedef List<edgeList> edgeListList;
}

#endif
