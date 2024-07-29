#ifndef uintSB_H
#define uintSB_H

#include "uint32.H"
#include "uint64.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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


MAXMIN(uint8_t, uint8_t, uint8_t)
MAXMIN(uint16_t, uint16_t, uint16_t)

MAXMIN(uint32_t, uint32_t, uint32_t)
MAXMIN(uint64_t, uint64_t, uint32_t)
MAXMIN(uint64_t, uint32_t, uint64_t)
MAXMIN(uint64_t, uint64_t, uint64_t)


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

unsigned int readUint(Istream&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
