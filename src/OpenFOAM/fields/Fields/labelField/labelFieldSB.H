// THE FIELD TEMPLATE IS NOW DEFINED FOR LABEL, AND BASIC ALGEBRAIC 
// MANIPULATIONS ARE ALSO DEFINED HERE.

#ifndef labelFieldSB_H
#define labelFieldSB_H

#include "Field.H"

#define TEMPLATE
#include "FieldFunctionsMSB.H"

namespace Foam
{

typedef Field<label> labelField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(label, label, mag)
BINARY_TYPE_OPERATOR(label, label, label, +, add)
BINARY_TYPE_OPERATOR(label, label, label, -, subtract)


template<>
tmp<labelField> labelField::component(const direction) const;

template<>
void component
(
    labelField& lf,
    const labelUList& f,
    const direction
);

template<>
void labelField::replace(const direction, const labelUList& lf);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "undefFieldFunctionsM.H"

#endif
