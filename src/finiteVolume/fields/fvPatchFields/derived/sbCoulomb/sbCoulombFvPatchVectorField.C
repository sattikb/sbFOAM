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
    alpha_(),
    motType_()
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
    alpha_(dict.lookup<scalar>("alpha")),
    motType_(dict.lookup<scalar>("motionType"))
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
    alpha_(ptf.alpha_),
    motType_(ptf.motType_)
{}


Foam::sbCoulombFvPatchVectorField::sbCoulombFvPatchVectorField
(
    const sbCoulombFvPatchVectorField& sbCoul,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(sbCoul, iF),
    omega_(sbCoul.omega_),
    alpha_(sbCoul.alpha_),
    motType_(sbCoul.motType_)
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

	tensor devTensor = muf*(gradUf + gradUf.T());

	vector stressWS = devTensor & nf;
	scalar tauW = mag(stressWS);

	scalar tauLim = pf*alpha_;
	
	    
	if (tauW>tauLim)
	{
       	    localSlipFaces++;

 	    scalar ry = p.deltaCoeffs()[faceI];
 	    scalar term1S = tauLim/(muf*ry);
 	    vector Utan = Uc - (Uc&nf)*nf;
	    vector tHat = Utan/max(mag(Utan),SMALL);

	    if(mag(Utan)<SMALL) Info<<"WHAT IS THIS"<<endl;
 	    Up[faceI] = Utan - (term1S * tHat);
	}
	else
	{
	    if(motType_==0)	//rotational
	    {
	    	const vector plateVel = omega_ ^ patch().Cf()[faceI];
	    	Up[faceI] = plateVel;
	    }
	    else		//translational
	    {
		Up[faceI] = omega_;
	    }
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
    writeEntry(os, "motionType", motType_);
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
