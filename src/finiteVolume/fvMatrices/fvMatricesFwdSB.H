// FWD DECLARATION OF ALL THE FvMatrix INSTANCES BASED OF FIELDTYPE

#ifndef fvMatricesFwdSB_H
#define fvMatricesFwdSB_H

#include "fieldTypes.H"	//CONTAINS THE FIELD TYPES SCALAR, VECTOR, TENSOR, ETC

namespace Foam
{

template<class Type>
class fvMatrix;

typedef fvMatrix<scalar> fvScalarMatrix;
typedef fvMatrix<vector> fvVectorMatrix;
typedef fvMatrix<sphericalTensor> fvSphericalTensorMatrix;
typedef fvMatrix<symmTensor> fvSymmTensorMatrix;
typedef fvMatrix<tensor> fvTensorMatrix;

} // End namespace Foam

#endif
