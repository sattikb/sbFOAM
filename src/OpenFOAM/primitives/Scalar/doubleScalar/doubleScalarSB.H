// FOR 32BIT (1/11/52) FLOATING POINT SCALARS

#ifndef doubleScalarSB_H
#define doubleScalarSB_H

#include "intSB.H"	// for the MAXMIN function used here.
#include "directionSB.H"

#include <limits>
#include <cmath>
using std::numeric_limits;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef double doubleScalar;

// Largest and smallest scalar values allowed in certain parts of the code.
static const doubleScalar doubleScalarVGreat
    = numeric_limits<doubleScalar>::max()/10;
static const doubleScalar doubleScalarVSmall
    = numeric_limits<doubleScalar>::min();

static const doubleScalar doubleScalarSmall
    = numeric_limits<doubleScalar>::epsilon();
static const doubleScalar doubleScalarGreat
    = 1.0/doubleScalarSmall;

static const doubleScalar doubleScalarRootVGreat = ::sqrt(doubleScalarVGreat);
static const doubleScalar doubleScalarRootVSmall = ::sqrt(doubleScalarVSmall);

static const doubleScalar doubleScalarRootGreat = ::sqrt(doubleScalarGreat);
static const doubleScalar doubleScalarRootSmall = ::sqrt(doubleScalarSmall);

static const doubleScalar doubleScalarNaN
    = numeric_limits<doubleScalar>::signaling_NaN();

//- Read whole of buf as a scalar. Return true if successful.
inline bool readScalar(const char* buf, doubleScalar& s)
{
    char* endPtr;
    s = strtod(buf, &endPtr);

    return (*endPtr == '\0');
}

#define Scalar doubleScalar
#define ScalarVGreat doubleScalarVGreat
#define ScalarVSmall doubleScalarVSmall
#define ScalarRootVGreat doubleScalarRootVGreat
#define ScalarRootVSmall doubleScalarRootVSmall
#define readScalar readDoubleScalar

inline Scalar mag(const Scalar s)
{
    return ::fabs(s);
}


#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return ::pow(Scalar(s), Scalar(e));           \
}

MAXMINPOW(Scalar, Scalar, Scalar)
MAXMINPOW(Scalar, Scalar, int)
MAXMINPOW(Scalar, int, Scalar)
MAXMINPOW(Scalar, Scalar, long)
MAXMINPOW(Scalar, long, Scalar)
MAXMINPOW(Scalar, Scalar, float)
MAXMINPOW(Scalar, float, Scalar)

#undef MAXMINPOW

#define transFunc(func)            \
inline Scalar func(const Scalar s) \
{                                  \
    return ::func(s);              \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline Scalar hypot(const Scalar x, const Scalar y)
{
    return ::hypot(x, y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    return ::jn(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    return ::yn(n, s);
}

#undef Scalar
#undef ScalarVGreat
#undef ScalarVSmall
#undef ScalarRootVGreat
#undef ScalarRootVSmall
#undef readScalar
#undef transFunc

} // End namespace Foam

#endif
