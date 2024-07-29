// DEPENDING ON WM_SP, WM_DP OR WM_LP IN BASHRC, ONE WILL BE CHOSEN AND 
// A SIMPLER NAME WILL BE ASSIGNED TO IT

#ifndef scalar_H
#define scalar_H

#include "floatScalarSB.H"	// FOR 32 BIT FLOATS ( 4 BYTES)
#include "doubleScalarSB.H"	// FOR 64 BIT FLOATS ( 8 BYTES)
#include "longDoubleScalarSB.H"	// FOR 80 BIT FLOATS (10 BYTES)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP)

// Define scalar as a float

namespace Foam
{
    typedef floatScalar scalar;

    static const scalar great = floatScalarGreat;
    static const scalar rootGreat = floatScalarRootGreat;
    static const scalar vGreat = floatScalarVGreat;
    static const scalar rootVGreat = floatScalarRootVGreat;
    static const scalar small = floatScalarSmall;
    static const scalar rootSmall = floatScalarRootSmall;
    static const scalar vSmall = floatScalarVSmall;
    static const scalar rootVSmall = floatScalarRootVSmall;
    static const scalar NaN = floatScalarNaN;
}

#elif defined(WM_DP)

// Define scalar as a double

namespace Foam
{
    typedef doubleScalar scalar;

    static const scalar great = doubleScalarGreat;
    static const scalar rootGreat = doubleScalarRootGreat;
    static const scalar vGreat = doubleScalarVGreat;
    static const scalar rootVGreat = doubleScalarRootVGreat;
    static const scalar small = doubleScalarSmall;
    static const scalar rootSmall = doubleScalarRootSmall;
    static const scalar vSmall = doubleScalarVSmall;
    static const scalar rootVSmall = doubleScalarRootVSmall;
    static const scalar NaN = doubleScalarNaN;
}

#elif defined(WM_LP)

// Define scalar as a long double

namespace Foam
{
    typedef longDoubleScalar scalar;

    static const scalar great = longDoubleScalarGreat;
    static const scalar rootGreat = longDoubleScalarRootGreat;
    static const scalar vGreat = longDoubleScalarVGreat;
    static const scalar rootVGreat = longDoubleScalarRootVGreat;
    static const scalar small = longDoubleScalarSmall;
    static const scalar rootSmall = longDoubleScalarRootSmall;
    static const scalar vSmall = longDoubleScalarVSmall;
    static const scalar rootVSmall = longDoubleScalarRootVSmall;
    static const scalar NaN = longDoubleScalarNaN;
}

#else

    #error "Precision not set, please set either WM_SP, WM_DP or WM_LP"

#endif

//- Deprecated limit constant for backward-compatibility
namespace Foam
{
    static const scalar GREAT = great;
    static const scalar ROOTGREAT = rootGreat;
    static const scalar VGREAT = vGreat;
    static const scalar ROOTVGREAT = rootVGreat;
    static const scalar SMALL = small;
    static const scalar ROOTSMALL = rootSmall;
    static const scalar VSMALL = vSmall;
    static const scalar ROOTVSMALL = rootVSmall;
}

// Additional global and transcendental functions

namespace Foam
{
    scalar readScalar(Istream& is);

    //- Inverse normalised incomplete gamma function
    scalar invIncGamma(const scalar a, const scalar P);

    //- Normalised upper incomplete gamma function
    scalar incGammaRatio_Q(const scalar a, const scalar x);

    //- Normalised lower incomplete gamma function
    scalar incGammaRatio_P(const scalar a, const scalar x);

    //- Upper incomplete gamma function
    scalar incGamma_Q(const scalar a, const scalar x);

    //- Lower incomplete gamma function
    scalar incGamma_P(const scalar a, const scalar x);
}

#endif
