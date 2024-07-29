// FORWARD DECLARATION OF SCALAR, VECTOR AND TENSOR FIELDS

#ifndef primitiveFieldsFwdSB_H
#define primitiveFieldsFwdSB_H

#include "fieldTypes.H"

namespace Foam
{

template<class Type> class Field;

typedef Field<label> labelField;
typedef Field<scalar> scalarField;
typedef Field<vector> vectorField;
typedef Field<sphericalTensor> sphericalTensorField;
typedef Field<symmTensor> symmTensorField;
typedef Field<tensor> tensorField;


} // End namespace Foam

#endif
