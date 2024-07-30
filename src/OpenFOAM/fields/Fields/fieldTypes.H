// INCLUDE THE HEADER FILES FOR ALL DIFFERENT KINDS OF FIELDS

#ifndef fieldTypes_H
#define fieldTypes_H

#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "sphericalTensorSB.H"
#include "symmTensor.H"
#include "tensor.H"
#include "triad.H"
#include "macros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define CAPITALIZE_scalar Scalar
#define CAPITALIZE_vector Vector
#define CAPITALIZE_sphericalTensor SphericalTensor
#define CAPITALIZE_symmTensor SymmTensor
#define CAPITALIZE_tensor Tensor

#define FOR_ALL_FIELD_TYPES(Macro, ...)                                        \
    Macro(scalar, __VA_ARGS__)                                                 \
    Macro(vector, __VA_ARGS__)                                                 \
    Macro(sphericalTensor, __VA_ARGS__)                                        \
    Macro(symmTensor, __VA_ARGS__)                                             \
    Macro(tensor, __VA_ARGS__)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
