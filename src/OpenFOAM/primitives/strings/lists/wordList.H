//DEFINES A LIST (1-d VECTOR) OF WORDS

#ifndef wordList_H
#define wordList_H

#include "word.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef UList<word> wordUList;

    typedef List<word> wordList;
}

#endif
