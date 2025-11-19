#ifndef sbPartialSlipMooneyFvPatchVectorField_H
#define sbPartialSlipMooneyFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sbPartialSlipMooneyFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sbPartialSlipMooneyFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
    vector omega_; 
    scalar alpha_;
    scalar m_;    

public:
    //- Runtime type information
    TypeName("sbPartialSlipMooney");

    // Constructors
    //- Construct from patch and internal field
    sbPartialSlipMooneyFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    sbPartialSlipMooneyFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given sbPartialSlipMooneyFvPatchVectorField onto a new patch
    sbPartialSlipMooneyFvPatchVectorField
    (
        const sbPartialSlipMooneyFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    sbPartialSlipMooneyFvPatchVectorField
    (
        const sbPartialSlipMooneyFvPatchVectorField&
    ) = delete;

    //- Copy constructor setting internal field reference
    sbPartialSlipMooneyFvPatchVectorField
    (
        const sbPartialSlipMooneyFvPatchVectorField&,
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
            new sbPartialSlipMooneyFvPatchVectorField(*this, iF)
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
