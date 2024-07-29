#include "longDoubleScalarSB.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define Scalar longDoubleScalar
#define ScalarVGreat longDoubleScalarVGreat
#define ScalarVSmall longDoubleScalarVSmall
#define ScalarRootVGreat longDoubleScalarRootVGreat
#define ScalarRootVSmall longDoubleScalarRootVSmall
#define readScalar readLongDoubleScalar
#include "Scalar.C"
#undef Scalar
#undef ScalarVGreat
#undef ScalarVSmall
#undef ScalarRootVGreat
#undef ScalarRootVSmall
#undef readScalar

// ************************************************************************* //
