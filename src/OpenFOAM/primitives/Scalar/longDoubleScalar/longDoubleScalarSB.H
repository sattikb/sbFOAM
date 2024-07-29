// FOR 80BIT (1/15/63) FLOATING POINT SCALARS

#ifndef longDoubleScalarSB_H
#define longDoubleScalarSB_H

#include "intSB.H"	 // for the MAXMIN function used here.
#include "directionSB.H"

#include <limits>
#include <cmath>
using std::numeric_limits;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef long double longDoubleScalar;

// Largest and smallest scalar values allowed in certain parts of the code.
static const longDoubleScalar longDoubleScalarVGreat
    = numeric_limits<longDoubleScalar>::max()/10;
static const longDoubleScalar longDoubleScalarVSmall
    = numeric_limits<longDoubleScalar>::min();

static const longDoubleScalar longDoubleScalarSmall
    = 1e3*numeric_limits<longDoubleScalar>::epsilon();
static const longDoubleScalar longDoubleScalarGreat
    = 1.0/longDoubleScalarSmall;

static const longDoubleScalar longDoubleScalarRootVGreat
    = ::sqrtl(longDoubleScalarVGreat);
static const longDoubleScalar longDoubleScalarRootVSmall
    = ::sqrtl(longDoubleScalarVSmall);

static const longDoubleScalar longDoubleScalarRootGreat
    = ::sqrtl(longDoubleScalarGreat);
static const longDoubleScalar longDoubleScalarRootSmall
    = ::sqrtl(longDoubleScalarSmall);

static const longDoubleScalar longDoubleScalarNaN
    = numeric_limits<longDoubleScalar>::signaling_NaN();

//- Read whole of buf as a scalar. Return true if successful.
inline bool readScalar(const char* buf, longDoubleScalar& s)
{
    char* endPtr;
    s = strtold(buf, &endPtr);

    return (*endPtr == '\0');
}

#define Scalar longDoubleScalar
#define ScalarVGreat longDoubleScalarVGreat
#define ScalarVSmall longDoubleScalarVSmall
#define ScalarRootVGreat longDoubleScalarRootVGreat
#define ScalarRootVSmall longDoubleScalarRootVSmall
#define readScalar readLongDoubleScalar

inline Scalar mag(const Scalar s)
{
    return ::fabsl(s);
}

#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline double pow(const type1 s, const type2 e)   \
{                                                 \
    return ::powl(Scalar(s), Scalar(e));          \
}

MAXMINPOW(Scalar, Scalar, Scalar)
MAXMINPOW(Scalar, Scalar, int)
MAXMINPOW(Scalar, int, Scalar)
MAXMINPOW(Scalar, Scalar, long)
MAXMINPOW(Scalar, long, Scalar)
MAXMINPOW(Scalar, Scalar, float)
MAXMINPOW(Scalar, float, Scalar)
MAXMINPOW(Scalar, Scalar, double)
MAXMINPOW(Scalar, double, Scalar)

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Scalar hypot(const Scalar x, const Scalar y)
{
    return ::hypotl(x, y);
}

inline Scalar atan2(const Scalar y, const Scalar x)
{
    return ::atan2l(y, x);
}

inline Scalar jn(const int n, const Scalar s)
{
    return ::jnl(n, s);
}

inline Scalar yn(const int n, const Scalar s)
{
    return ::ynl(n, s);
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
