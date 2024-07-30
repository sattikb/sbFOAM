// A SYMMETRIC TENSOR CREATED OUT OF SCALARS. ONLY NEED INFORMATION ABOUT 6
// ELEMENTS FOR A 3X3 MATRIX

#ifndef symmTensorSB_H
#define symmTensorSB_H

#include "SymmTensorSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef SymmTensor<scalar> symmTensor;

//- Data associated with symmTensor type are contiguous
template<>
inline bool contiguous<symmTensor>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
