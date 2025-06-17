#include "opeSlipSBFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::opeSlipSBFvPatchVectorField::opeSlipSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    k1_(),
    k2_(),
    C_(),
    vp_(),
    sigma_()
{}


Foam::opeSlipSBFvPatchVectorField::opeSlipSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    k1_(dict.lookup<scalar>("k1")),
    k2_(dict.lookup<scalar>("k2")),
    C_(dict.lookup<scalar>("C")),
    vp_(dict.lookup<vector>("surfaceVelocity")),
    sigma_(dict.lookup<vector>("sigma"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::opeSlipSBFvPatchVectorField::opeSlipSBFvPatchVectorField
(
    const opeSlipSBFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    k1_(ptf.k1_),
    k2_(ptf.k2_),
    C_(ptf.C_),
    vp_(ptf.vp_),
    sigma_(ptf.sigma_)
{}


Foam::opeSlipSBFvPatchVectorField::opeSlipSBFvPatchVectorField
(
    const opeSlipSBFvPatchVectorField& oswSB,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(oswSB, iF),
    k1_(oswSB.k1_),
    k2_(oswSB.k2_),
    C_(oswSB.C_),
    vp_(oswSB.vp_),
    sigma_(oswSB.sigma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::opeSlipSBFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the velocity field from the database
//    const volVectorField& U = db().lookupObject<volVectorField>("U");
//    const vectorField& Uc = U.boundaryField()[patch().index()];
//    const vectorField n = patch().nf(); // Normal vectors pointing out of the domain

    const volScalarField& pField = db().lookupObject<volScalarField>("p");
    const fvPatchScalarField& pp = pField.boundaryField()[patch().index()];
    // Initialize patch velocity field
    vectorField& Up = *this;

    const scalar sigmaG  = sigma_[0];
    const scalar sigmaNS = sigma_[1];
    const scalar vMax    = 10;

    forAll(Up, faceI)
    {
	const scalar sigmaStar = (pp[faceI] - sigmaNS) / (sigmaG - sigmaNS);
        const scalar vc        = C_ * sigmaStar;
	const scalar vpStar    = mag(vp_) / vMax;
	const vector dirVec   = (vpStar > SMALL) ? vp_ / vpStar : vector::zero;

	const scalar term = -k2_ * vc;
        scalar xi = 0.0;
        if (term<15) xi  = 1.0 / (1.0 + exp(term));

        const scalar num   = 1.0 + exp(k1_ * (vc - vpStar));
        const scalar den   = 1.0 + exp(k2_ * vc);
        const scalar frac  = log(num) / (log(den) + SMALL);

        scalar vw = vc * xi * (1.0 - frac);
        if (vw > vc) vw = 0.0;

        Up[faceI] = vw * dirVec;

    }
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::opeSlipSBFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "k1", k1_);
    writeEntry(os, "k2", k2_);
    writeEntry(os, "C", C_);
    writeEntry(os, "surfaceVelocity", vp_);
    writeEntry(os, "sigma", sigma_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        opeSlipSBFvPatchVectorField
    );
}
