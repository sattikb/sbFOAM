#include "sbCoulombFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sbCoulombFvPatchVectorField::sbCoulombFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    omega_(),
    alpha_()
{}


Foam::sbCoulombFvPatchVectorField::sbCoulombFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    omega_(dict.lookup<vector>("shaftOmega")),
    alpha_(dict.lookup<scalar>("alpha"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::sbCoulombFvPatchVectorField::sbCoulombFvPatchVectorField
(
    const sbCoulombFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_),
    alpha_(ptf.alpha_)
{}


Foam::sbCoulombFvPatchVectorField::sbCoulombFvPatchVectorField
(
    const sbCoulombFvPatchVectorField& sbCoul,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(sbCoul, iF),
    omega_(sbCoul.omega_),
    alpha_(sbCoul.alpha_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sbCoulombFvPatchVectorField::updateCoeffs()
{
    Info<<"SATTIK IN COULOMB TYPE SLIP MODEL."<<endl;
    label localProcessedFaces = 0;
    label localSlipFaces = 0;

    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();
    const vectorField n = p.nf();
    const vectorField Sfp = p.Sf();

    const volVectorField& U = db().lookupObject<volVectorField>("U");
    const volScalarField& mu = db().lookupObject<volScalarField>("thermo:mu");
    const volScalarField& pField = db().lookupObject<volScalarField>("p");

    const fvPatchScalarField& pp = pField.boundaryField()[p.index()];
    const fvPatchScalarField& mup = mu.boundaryField()[p.index()];

    tensorField gradU = fvc::grad(U)().boundaryField()[p.index()];

    vectorField& Up = *this;
    forAll(Up, faceI)
    {
	localProcessedFaces++;
        scalar muf = mup[faceI];
        scalar pf = pp[faceI];

 	label own = p.faceCells()[faceI];
 	vector Uc = U[own];
 	vector Uf = U[faceI];

        const vector& nf = n[faceI];
        const tensor& gradUf = gradU[faceI];

	tensor tauf = muf*(gradUf + gradUf.T()) - (2.0/3.0)*muf*(tr(gradUf))*I;
	tensor sigmaf = tauf - pf*I;

// 1ST METHOD TO COMPUTE WALL SHEAR STRESS
	vector traction = sigmaf & nf;
	scalar stressWN = traction & nf;
	vector stressWS = traction  - stressWN*nf;
	scalar tauW = mag(stressWS);
	vector tauHat =  stressWS / max(tauW,SMALL);	

// 2ND METHOD TO COMPUTE WALL SHEAR STRESS
	vector stressWS2 = nf & tauf;
	scalar tauW2 = mag(stressWS2);

// 3RD METHOD TO COMPUTE WALL SHEAR STRESS
	vector UcTan = Uc - (Uc&nf)*nf;
 	scalar dsInv = p.deltaCoeffs()[faceI];
	vector stressWS3 = muf*(-UcTan)*dsInv;
	scalar tauW3 = mag(stressWS3);
	//vector tauHat =  stressWS3 / max(tauW3,SMALL);	

	scalar tauLim = pf*alpha_;
	
	    
	if (tauW>tauLim)
	{
       	    localSlipFaces++;

 	    scalar ry = p.deltaCoeffs()[faceI];
 	    scalar term1 = tauLim/(muf*ry);
 	    scalar term2 = Uc & tauHat;
 	    Up[faceI] = (term1 + term2) * tauHat;
 	   // Up[faceI] = (term1) * tauHat + UcTan;
	}
	else
	{
	    const vector plateVel = omega_ ^ patch().Cf()[faceI];
	    Up[faceI] = plateVel;
	}

    }
    
    label globalProcessed = returnReduce(localProcessedFaces, sumOp<label>());
    label globalSlip      = returnReduce(localSlipFaces,      sumOp<label>());
    Info << "=== COULOMB SLIP BC ACTIVE === "
         << "processed " << globalProcessed
         << " faces globally, slip on " << globalSlip
	 << " faces this update (patch: " << patch().name() << ")"
	 << endl;

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::sbCoulombFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "shaftOmega", omega_);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sbCoulombFvPatchVectorField
    );
}
