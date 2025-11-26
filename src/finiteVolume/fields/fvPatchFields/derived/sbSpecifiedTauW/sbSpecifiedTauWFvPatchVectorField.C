#include "sbSpecifiedTauWFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

// Constructors
sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    f_()
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    f_(dict.lookup<scalar>("f"))
{
    // Initialize with zero gradient, will be updated in updateCoeffs
    gradient() = vectorField(p.size(), vector::zero);
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }
}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const sbSpecifiedTauWFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
    f_(ptf.f_)
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const sbSpecifiedTauWFvPatchVectorField& sbstw,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(sbstw, iF),
    f_(sbstw.f_)
{}

// Member Functions
void sbSpecifiedTauWFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get pressure field
    const volScalarField& p = db().lookupObject<volScalarField>("p");
    const fvPatchScalarField& pp = p.boundaryField()[patch().index()];

    // Compute gradient: du/dn = f * p
    gradient() = f_ * pp * vector::zero;

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void sbSpecifiedTauWFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    writeEntry(os, "f", f_);
}

makePatchTypeField
(
    fvPatchVectorField,
    sbSpecifiedTauWFvPatchVectorField
);

} // End namespace Foam
