#ifndef sbSpecifiedTauWFvPatchVectorField_H
#define sbSpecifiedTauWFvPatchVectorField_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class sbSpecifiedTauWFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
private:
    // Private Data
    scalar f_;  // Coefficient vector for du/dn = f * p

public:
    //- Runtime type information
    TypeName("sbSpecifiedTauW");

    // Constructors
    sbSpecifiedTauWFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    sbSpecifiedTauWFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping
    sbSpecifiedTauWFvPatchVectorField
    (
        const sbSpecifiedTauWFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy setting internal field reference
    sbSpecifiedTauWFvPatchVectorField
    (
        const sbSpecifiedTauWFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new sbSpecifiedTauWFvPatchVectorField(*this, iF)
        );
    }

    // Member Functions
    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
