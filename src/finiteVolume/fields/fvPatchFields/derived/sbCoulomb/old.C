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

    // Get the velocity field from the database
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

        const vector& nf = n[faceI];
        const tensor& gradUf = gradU[faceI];

	tensor tauf = muf*(gradUf + gradUf.T()) - (2.0/3.0)*muf*(tr(gradUf))*I;
	tensor sigmaf = tauf - pf*I;

	vector traction = sigmaf & nf;
	scalar stressWN = traction & nf;
	vector stressWS = traction  - stressWN*nf;
	scalar tauW = mag(stressWS);
	vector tauHat =  stressWS / max(tauW,SMALL);	

	const vector plateVel = omega_ ^ patch().Cf()[faceI];
	scalar tauLim = pf*alpha_;


	Up[faceI] = plateVel;
	if (tauW>tauLim)
	{
       	    localSlipFaces++;

 	    label own = p.faceCells()[faceI];
 	    vector Uc = U[own];
 	    const scalarField& ry = p.deltaCoeffs();
// 	    vector fromESI = tauHat*(tauHat & (stressWS/(muf*ry[faceI]) + Uc));
 	    vector newStress = tauHat * tauLim;
 	    vector fromESI = newStress/(muf*ry[faceI]) + (Uc & tauHat)*tauHat;
	    //Info<<"First term: "<<newStress/(muf*ry[faceI])<<" Second term: "<<(Uc & tauHat)*tauHat<<endl;
 
 	    
 	    Up[faceI] = fromESI;
	}
//	Up[faceI] = plateVel;
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
