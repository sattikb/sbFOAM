#include "opeRotSlipSBFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::opeRotSlipSBFvPatchVectorField::opeRotSlipSBFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    k1_(),
    k2_(),
    C_(),
    omegap_(),
    sigma_()
{}


Foam::opeRotSlipSBFvPatchVectorField::opeRotSlipSBFvPatchVectorField
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
    omegap_(dict.lookup<vector>("surfaceOmega")),
    sigma_(dict.lookup<vector>("sigma"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::opeRotSlipSBFvPatchVectorField::opeRotSlipSBFvPatchVectorField
(
    const opeRotSlipSBFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    k1_(ptf.k1_),
    k2_(ptf.k2_),
    C_(ptf.C_),
    omegap_(ptf.omegap_),
    sigma_(ptf.sigma_)
{}


Foam::opeRotSlipSBFvPatchVectorField::opeRotSlipSBFvPatchVectorField
(
    const opeRotSlipSBFvPatchVectorField& orswSB,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(orswSB, iF),
    k1_(orswSB.k1_),
    k2_(orswSB.k2_),
    C_(orswSB.C_),
    omegap_(orswSB.omegap_),
    sigma_(orswSB.sigma_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::opeRotSlipSBFvPatchVectorField::updateCoeffs()
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
    	const scalar vc = C_ * sigmaStar;
        const vector plateVel = omegap_ ^ patch().Cf()[faceI];
	const scalar vp       = mag(plateVel) / vMax; 
	const vector dirVec   = (vp > SMALL) ? plateVel / vp : vector::zero;
    	
        const scalar xi  = 1.0 / (1.0 + exp(-k2_ * vc));
//	Info<<"SATTIK p, vc, xi  VALUES: "<<pp[faceI]<<" :: "<<vc<<" :: "<<xi<<endl;

        const scalar num   = 1.0 + exp(k1_ * (vc - vp));
        const scalar den   = 1.0 + exp(k2_ * vc);
        const scalar frac  = log(num) / (log(den) + SMALL);
//	Info<<"SATTIK num,den  VALUES: "<<num<<" :: "<<den<<endl;
        
	scalar vw = vc * xi * (1.0 - frac);
        if (vw > vc) vw = 0.0;
//	Info<<"SATTIK vw"<<vw<<endl;

//	Info<<"SATTIK pressure, vcScalar  VALUES: "<<pp[faceI]<<" : "<<vcCCC<<endl;
//	Info<<"SATTIK vw,vp  VALUES: "<<vw<<"::"<<vp<<endl;
    	Up[faceI] = vw * dirVec;
    }
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::opeRotSlipSBFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "k1", k1_);
    writeEntry(os, "k2", k2_);
    writeEntry(os, "C", C_);
    writeEntry(os, "surfaceOmega", omegap_);
    writeEntry(os, "sigma", sigma_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        opeRotSlipSBFvPatchVectorField
    );
}
