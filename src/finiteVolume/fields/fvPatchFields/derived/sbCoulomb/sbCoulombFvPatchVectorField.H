#ifndef sbCoulombFvPatchVectorField_H
#define sbCoulombFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sbCoulombFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sbCoulombFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
    vector omega_; 
    scalar alpha_;

public:
    //- Runtime type information
    TypeName("sbCoulomb");

    // Constructors
    //- Construct from patch and internal field
    sbCoulombFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    sbCoulombFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given sbCoulombFvPatchVectorField onto a new patch
    sbCoulombFvPatchVectorField
    (
        const sbCoulombFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    sbCoulombFvPatchVectorField
    (
        const sbCoulombFvPatchVectorField&
    ) = delete;

    //- Copy constructor setting internal field reference
    sbCoulombFvPatchVectorField
    (
        const sbCoulombFvPatchVectorField&,
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
            new sbCoulombFvPatchVectorField(*this, iF)
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
