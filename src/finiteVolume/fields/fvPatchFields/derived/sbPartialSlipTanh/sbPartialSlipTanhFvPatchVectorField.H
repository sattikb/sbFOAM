#ifndef sbPartialSlipTanhFvPatchVectorField_H
#define sbPartialSlipTanhFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sbPartialSlipTanhFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sbPartialSlipTanhFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
    vector omega_; 
    scalar tauStar_;
    scalar C_; 
    scalar m_;    

public:
    //- Runtime type information
    TypeName("sbPartialSlipTanh");

    // Constructors
    //- Construct from patch and internal field
    sbPartialSlipTanhFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    sbPartialSlipTanhFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given sbPartialSlipTanhFvPatchVectorField onto a new patch
    sbPartialSlipTanhFvPatchVectorField
    (
        const sbPartialSlipTanhFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    sbPartialSlipTanhFvPatchVectorField
    (
        const sbPartialSlipTanhFvPatchVectorField&
    ) = delete;

    //- Copy constructor setting internal field reference
    sbPartialSlipTanhFvPatchVectorField
    (
        const sbPartialSlipTanhFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    //- Create a clone with internal field reference
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new sbPartialSlipTanhFvPatchVectorField(*this, iF)
        );
    }

    // Member Functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif
