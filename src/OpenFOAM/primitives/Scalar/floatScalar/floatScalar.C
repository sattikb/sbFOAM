#include "floatScalarSB.H"	//FOR 32BIT SCALARS.	
#include "IOstreams.H"

#include <sstream>


#define Scalar floatScalar
#define ScalarVGreat floatScalarVGreat
#define ScalarVSmall floatScalarVSmall
#define ScalarRootVGreat floatScalarRootVGreat
#define ScalarRootVSmall floatScalarRootVSmall
#define readScalar readFloatScalar
#include "Scalar.C"
#undef Scalar
#undef ScalarVSmall
#undef ScalarVSmall
#undef ScalarRootVGreat
#undef ScalarRootVSmall
#undef readScalar
