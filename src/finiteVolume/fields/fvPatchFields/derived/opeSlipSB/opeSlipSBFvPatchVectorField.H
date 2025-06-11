#ifndef opeSlipSBFvPatchVectorField_H
#define opeSlipSBFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class opeSlipSBFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class opeSlipSBFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
    scalar k1_;     // First constant (k1 = 50)
    scalar k2_;     // Second constant (k2 = 50)
    scalar C_;      // Constant (C = 1)
    vector vp_;     // Constant (C = 1)
    vector sigma_;  // User-specified sigma value

public:
    //- Runtime type information
    TypeName("opeSlipSB");

    // Constructors
    //- Construct from patch and internal field
    opeSlipSBFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    opeSlipSBFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given opeSlipSBFvPatchVectorField onto a new patch
    opeSlipSBFvPatchVectorField
    (
        const opeSlipSBFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    opeSlipSBFvPatchVectorField
    (
        const opeSlipSBFvPatchVectorField&
    ) = delete;

    //- Copy constructor setting internal field reference
    opeSlipSBFvPatchVectorField
    (
        const opeSlipSBFvPatchVectorField&,
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
            new opeSlipSBFvPatchVectorField(*this, iF)
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
