// THE FIELD TEMPLATE IS NOW SPECIALIZED FOR SPHERICAL TENSORS.

#ifndef sphericalTensorFieldSB_H
#define sphericalTensorFieldSB_H

#include "scalarField.H"
#include "sphericalTensorSB.H"

#define TEMPLATE
#include "FieldFunctionsMSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef Field<sphericalTensor> sphericalTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, sphericalTensor, tr)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, sph)
UNARY_FUNCTION(scalar, sphericalTensor, det)
UNARY_FUNCTION(sphericalTensor, sphericalTensor, inv)

BINARY_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)
BINARY_TYPE_OPERATOR(sphericalTensor, scalar, sphericalTensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "undefFieldFunctionsM.H"

#endif
