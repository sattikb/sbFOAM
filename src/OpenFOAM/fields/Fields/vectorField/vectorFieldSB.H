// THE FIELD TEMPLATE IS NOW TYPEDEFFED FOR VECTORS.

#ifndef vectorFieldSB_H
#define vectorFieldSB_H

#include "scalarField.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<vector> vectorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<vectorField> zip
(
    const tmp<scalarField>& x,
    const tmp<scalarField>& y,
    const tmp<scalarField>& z
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
