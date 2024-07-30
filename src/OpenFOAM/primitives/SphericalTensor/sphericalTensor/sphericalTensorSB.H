// THE SPHERICALTENSOR TEMPLATE IS USED FOR SCALARS.

#ifndef sphericalTensorSB_H
#define sphericalTensorSB_H

#include "IdentitySB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef SphericalTensor<scalar> sphericalTensor;

static const sphericalTensor oneThirdI(1.0/3.0);
static const sphericalTensor twoThirdsI(2.0/3.0);

//- Specify data associated with sphericalTensor type are contiguous
template<>
inline bool contiguous<sphericalTensor>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif
