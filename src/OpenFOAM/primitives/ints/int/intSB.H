//READS INTEGER, AND DEFINES MIN AND MAX FUNCTIONS FOR ALL INTEGER TYPES

#ifndef intSB_H
#define intSB_H

#include "int32.H"
#include "int64.H"

namespace Foam
{

#define MAXMIN(retType, type1, type2)              \
                                                   \
inline retType max(const type1 s1, const type2 s2) \
{                                                  \
    return (s1 > s2)? s1: s2;                      \
}                                                  \
                                                   \
inline retType min(const type1 s1, const type2 s2) \
{                                                  \
    return (s1 < s2)? s1: s2;                      \
}


MAXMIN(int8_t, int8_t, int8_t)
MAXMIN(int16_t, int16_t, int16_t)

MAXMIN(int32_t, int32_t, int32_t)
MAXMIN(int64_t, int64_t, int32_t)
MAXMIN(int64_t, int32_t, int64_t)
MAXMIN(int64_t, int64_t, int64_t)


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int readInt(Istream&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
