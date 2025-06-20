#include "fixedShearStressFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tau0_(Zero)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
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


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tau0_(ptf.tau0_)
{}


Foam::fixedShearStressFvPatchVectorField::fixedShearStressFvPatchVectorField
(
    const fixedShearStressFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    tau0_(ptf.tau0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedShearStressFvPatchVectorField::updateCoeffs()
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

    scalarField nuEff(turbModel.nuEff(patch().index()));

    const vectorField Uc(patchInternalField());

    vector tauHat = tau0_/(mag(tau0_) + rootVSmall);

    Info<<"SATTIK NUEFF: "<<nuEff<<endl;

    const scalarField& ry = patch().deltaCoeffs();

    operator==(tauHat*(tauHat & (tau0_*(1.0/(ry*nuEff)) + Uc)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::fixedShearStressFvPatchVectorField::write(Ostream& os) const
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
        fixedShearStressFvPatchVectorField
    );
}
