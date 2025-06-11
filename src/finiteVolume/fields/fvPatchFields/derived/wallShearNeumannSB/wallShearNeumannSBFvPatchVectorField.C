#include "wallShearNeumannSBFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "symmTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallShearNeumannSBFvPatchVectorField::wallShearNeumannSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), Zero),
    mu_(0.0)
{}


wallShearNeumannSBFvPatchVectorField::wallShearNeumannSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(dict.lookup("traction")),
    mu_(dict.lookup("mu"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    updateCoeffs();
}


wallShearNeumannSBFvPatchVectorField::wallShearNeumannSBFvPatchVectorField
(
    const wallShearNeumannSBFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
    traction_(mapper(ptf.traction_)),
    mu_(ptf.mu_)
{}


wallShearNeumannSBFvPatchVectorField::wallShearNeumannSBFvPatchVectorField
(
    const wallShearNeumannSBFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf),
    traction_(ptf.traction_),
    mu_(ptf.mu_)
{}


wallShearNeumannSBFvPatchVectorField::wallShearNeumannSBFvPatchVectorField
(
    const wallShearNeumannSBFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF),
    traction_(ptf.traction_),
    mu_(ptf.mu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wallShearNeumannSBFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get face unit normals
    const vectorField n(patch().nf());

    // Identity tensor
    const tensor I = tensor::I;

    // Compute gradient to satisfy τ · n = traction
    // τ = μ*(gradU + gradU.T) => gradU = (τ / (2μ)) · (I - n⊗n)
    gradient() = (traction_/(2.0*mu_)) & (I - sqr(n));

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void wallShearNeumannSBFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    os.writeEntry("mu", mu_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    wallShearNeumannSBFvPatchVectorField
);

} // End namespace Foam
