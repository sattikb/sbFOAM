#include "wallShearSBFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallShearSBFvPatchVectorField::wallShearSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tau0_(Zero)
{}


Foam::wallShearSBFvPatchVectorField::wallShearSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    tau0_(dict.lookupOrDefault<vector>("tau", Zero))
{
    fvPatchField<vector>::operator=(patchInternalField());
}


Foam::wallShearSBFvPatchVectorField::wallShearSBFvPatchVectorField
(
    const wallShearSBFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tau0_(ptf.tau0_)
{}


Foam::wallShearSBFvPatchVectorField::wallShearSBFvPatchVectorField
(
    const wallShearSBFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    tau0_(ptf.tau0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallShearSBFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    const volScalarField& rhoField = db().lookupObject<volScalarField>("rho");
    const scalarField rho(rhoField.boundaryField()[patch().index()]);
    
    scalarField nuEff(turbModel.nuEff(patch().index()));
    scalarField muEff = nuEff * rho;

    const vectorField Uc(patchInternalField());

    vector tauHat = tau0_/(mag(tau0_) + rootVSmall);

    const scalarField& ry = patch().deltaCoeffs();

    operator==(tauHat*(tauHat & (tau0_*(1.0/(ry*muEff)) + Uc)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::wallShearSBFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "tau", tau0_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        wallShearSBFvPatchVectorField
    );
}
