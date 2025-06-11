#ifndef opeRotSlipSBFvPatchVectorField_H
#define opeRotSlipSBFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class opeRotSlipSBFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class opeRotSlipSBFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
    scalar k1_;     // First constant (k1 = 50)
    scalar k2_;     // Second constant (k2 = 50)
    scalar C_;      // Constant (C = 1)
    vector omegap_;     // Constant (C = 1)
    vector sigma_;  // User-specified sigma value

public:
    //- Runtime type information
    TypeName("opeRotSlipSB");

    // Constructors
    //- Construct from patch and internal field
    opeRotSlipSBFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    opeRotSlipSBFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given opeRotSlipSBFvPatchVectorField onto a new patch
    opeRotSlipSBFvPatchVectorField
    (
        const opeRotSlipSBFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    opeRotSlipSBFvPatchVectorField
    (
        const opeRotSlipSBFvPatchVectorField&
    );

    //- Copy constructor setting internal field reference
    opeRotSlipSBFvPatchVectorField
    (
        const opeRotSlipSBFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    //- Create a clone
    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new opeRotSlipSBFvPatchVectorField(*this)
        );
    }

    //- Create a clone with internal field reference
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new opeRotSlipSBFvPatchVectorField(*this, iF)
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
