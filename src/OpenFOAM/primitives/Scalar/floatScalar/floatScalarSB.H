// FOR 32BIT (1/8/23) FLOATING POINT SCALARS

#ifndef floatScalarSB_H
#define floatScalarSB_H

#include "intSB.H"  // for the MAXMIN function used here. 
#include "directionSB.H"

#include <limits>
#include <cmath>
using std::numeric_limits;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef float floatScalar;

// Largest and smallest scalar values allowed in certain parts of the code
static const floatScalar floatScalarVGreat
    = numeric_limits<floatScalar>::max()/10;
static const floatScalar floatScalarVSmall
    = numeric_limits<floatScalar>::min();

static const floatScalar floatScalarSmall
    = numeric_limits<floatScalar>::epsilon();
static const floatScalar floatScalarGreat
    = 1.0/floatScalarSmall;

static const floatScalar floatScalarRootVGreat = ::sqrt(floatScalarVGreat);
static const floatScalar floatScalarRootVSmall = ::sqrt(floatScalarVSmall);

static const floatScalar floatScalarRootGreat = ::sqrt(floatScalarGreat);
static const floatScalar floatScalarRootSmall = ::sqrt(floatScalarSmall);

static const floatScalar floatScalarNaN
    = numeric_limits<floatScalar>::signaling_NaN();

//- Read whole of buf as a scalar. Return true if successful.
inline bool readScalar(const char* buf, floatScalar& s)
{
    char* endPtr;
    s = strtof(buf, &endPtr);

    return (*endPtr == '\0');
}

#define Scalar floatScalar
#define ScalarVGreat floatScalarVGreat
#define ScalarVSmall floatScalarVSmall
#define ScalarRootVGreat floatScalarRootVGreat
#define ScalarRootVSmall floatScalarRootVSmall
#define readScalar readFloatScalar

inline Scalar mag(const Scalar s)
{
    return ::fabsf(s);
}


#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return ::powf(s, e);                          \
}

MAXMINPOW(Scalar, Scalar, Scalar)
MAXMINPOW(Scalar, Scalar, int)
MAXMINPOW(Scalar, int, Scalar)
MAXMINPOW(Scalar, Scalar, long)
MAXMINPOW(Scalar, long, Scalar)

#undef MAXMINPOW


#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func##f(s);           \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Scalar hypot(const Scalar x, const Scalar y)
{
    return ::hypotf(x, y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2f(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    return ::jnf(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    return ::ynf(n, s);
}

#undef Scalar
#undef ScalarVGreat
#undef ScalarVSmall
#undef ScalarRootVGreat
#undef ScalarRootVSmall
#undef readScalar
#undef transFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
