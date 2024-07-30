// WILL BE USED TO CREATE A TENSOR MADE UP OF SCALARS

#ifndef tensor_H
#define tensor_H

#include "vector.H"
#include "sphericalTensorSB.H"
#include "symmTensorSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Tensor<scalar> tensor;

vector eigenValues(const tensor& T);
vector eigenVector
(
    const tensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
);
tensor eigenVectors(const tensor& T, const vector& lambdas);
tensor eigenVectors(const tensor& T);

vector eigenValues(const symmTensor& T);
vector eigenVector
(
    const symmTensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
);
tensor eigenVectors(const symmTensor& T, const vector& lambdas);
tensor eigenVectors(const symmTensor& T);

//- Data associated with tensor type are contiguous
template<>
inline bool contiguous<tensor>() {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
