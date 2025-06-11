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

    // Initialize patch velocity field
    vectorField& Up = *this;

    const vector vc = C_ * sigma_;

    forAll(Up, faceI)
    {
 

    	vector vw = vector::zero;
    	scalar nn = pTraits<vector>::nComponents;
    	for (direction ii = 0; ii< nn; ii++)
    	{
    	    const scalar vcComp  = vc.component(ii);
    	    const scalar vpComp  = vp_.component(ii);

    		const scalar xiComp  = 1.0 / (1.0 + exp(-k2_ * vcComp));

    	    const scalar num1   = 1.0 + exp(k1_ * (vcComp - vpComp));
    	    const scalar den1   = 1.0 + exp(k2_ * vcComp);
    	    const scalar frac1  = log(num1) / log(den1);
    	    scalar vwComp = vcComp * xiComp * (1.0 - frac1);

    	    if (vwComp > vcComp) vwComp = 0.0;
    	    vw.component(ii) = vwComp;
    	}

    	Up[faceI] = vw;
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
