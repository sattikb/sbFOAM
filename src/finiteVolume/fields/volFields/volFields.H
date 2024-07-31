#ifndef volFields_H
#define volFields_H

#include "GeometricFieldsSB.H"
#include "calculatedFvPatchFieldsSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<GeometricField<scalar, fvPatchField, volMesh>>
GeometricField<scalar, fvPatchField, volMesh>::component
(
    const direction
) const;

template<>
void GeometricField<scalar, fvPatchField, volMesh>::replace
(
    const direction,
    const GeometricField<scalar, fvPatchField, volMesh>& sf
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "volFieldsI.H"

#endif
