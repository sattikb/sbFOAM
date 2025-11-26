#ifndef sbPartialSlipOnlyTauFvPatchVectorField_H
#define sbPartialSlipOnlyTauFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "movingWallVelocityFvPatchVectorField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sbPartialSlipOnlyTauFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sbPartialSlipOnlyTauFvPatchVectorField
:
    public movingWallVelocityFvPatchVectorField // fixedValueFvPatchVectorField
{
private:
    // Private Data
    vector omega_; 
    scalar alpha_; 
    scalar m_;    

public:
    //- Runtime type information
    TypeName("sbPartialSlipOnlyTau");

    // Constructors
    //- Construct from patch and internal field
    sbPartialSlipOnlyTauFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    sbPartialSlipOnlyTauFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given sbPartialSlipOnlyTauFvPatchVectorField onto a new patch
    sbPartialSlipOnlyTauFvPatchVectorField
    (
        const sbPartialSlipOnlyTauFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    sbPartialSlipOnlyTauFvPatchVectorField
    (
        const sbPartialSlipOnlyTauFvPatchVectorField&
    ) = delete;

    //- Copy constructor setting internal field reference
    sbPartialSlipOnlyTauFvPatchVectorField
    (
        const sbPartialSlipOnlyTauFvPatchVectorField&,
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
            new sbPartialSlipOnlyTauFvPatchVectorField(*this, iF)
        );
    }

    // Member Functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;

    virtual bool moving() const; 
};


} // End namespace Foam

#endif
