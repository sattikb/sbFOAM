// FORWARD DECLARATION OF THE VARIOUS FIELD TYPES.

#ifndef volFieldsFwd_H
#define volFieldsFwd_H

#include "fieldTypes.H"

namespace Foam
{

class volMesh;

template<class Type>
class fvPatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

template<class Type>
using VolField = GeometricField<Type, fvPatchField, volMesh>;

typedef GeometricField<scalar,          fvPatchField, volMesh> volScalarField;
typedef GeometricField<vector,          fvPatchField, volMesh> volVectorField;
typedef GeometricField<sphericalTensor, fvPatchField, volMesh> volSphericalTensorField;
typedef GeometricField<symmTensor,      fvPatchField, volMesh> volSymmTensorField;
typedef GeometricField<tensor,          fvPatchField, volMesh> volTensorField;

} // End namespace Foam

#endif
