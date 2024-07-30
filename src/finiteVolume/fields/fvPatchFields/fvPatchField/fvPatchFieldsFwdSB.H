// THE TEMPLATED FVPATCHFIELD IS NOW TYPEDEFFED WITH VARIOUS DATA TYPES.

#ifndef fvPatchFieldsFwdSB_H
#define fvPatchFieldsFwdSB_H

//#include "fieldTypes.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type> class fvPatchField;

typedef fvPatchField<scalar> fvPatchScalarField;
typedef fvPatchField<vector> fvPatchVectorField;
typedef fvPatchField<sphericalTensor> fvPatchSphericalTensorField;
typedef fvPatchField<symmTensor> fvPatchSymmTensorField;
typedef fvPatchField<tensor> fvPatchTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
