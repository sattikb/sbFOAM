#include "sbPartialSlipMooneyFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sbPartialSlipMooneyFvPatchVectorField::sbPartialSlipMooneyFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    omega_(),
    alpha_(),
    m_()
{}


Foam::sbPartialSlipMooneyFvPatchVectorField::sbPartialSlipMooneyFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    omega_(dict.lookup<vector>("shaftOmega")),
    alpha_(dict.lookup<scalar>("alpha")),
    m_(dict.lookup<scalar>("m"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::sbPartialSlipMooneyFvPatchVectorField::sbPartialSlipMooneyFvPatchVectorField
(
    const sbPartialSlipMooneyFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_),
    alpha_(ptf.alpha_),
    m_(ptf.m_)
{}


Foam::sbPartialSlipMooneyFvPatchVectorField::sbPartialSlipMooneyFvPatchVectorField
(
    const sbPartialSlipMooneyFvPatchVectorField& pstauSB,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pstauSB, iF),
    omega_(pstauSB.omega_),
    alpha_(pstauSB.alpha_),
    m_(pstauSB.m_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sbPartialSlipMooneyFvPatchVectorField::updateCoeffs()
{
    Info<<"SATTIK IN MOONEY PARTIAL SLIP."<<endl;
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
	scalar tauN = - tauW / stressWN;	

	const vector plateVel = omega_ ^ patch().Cf()[faceI];


	Up[faceI] = plateVel;
	if (tauN>0)
	{
		vector dir = -stressWS/tauW;
	
//		Info<<"------ value of dir is: "<<dir<<endl;
		scalar phi = (alpha_)/stressWN*pow(tauN, m_);
//		Info<<"------ value of tanNorm and tanValFrac is: "<<tauN<<" and "<<tanVal<<endl;
		
		label own = p.faceCells()[faceI];
		vector cellVel = U[own];
		scalar magCellVel = mag(cellVel);
		
		Up[faceI] = (1-phi)*plateVel + (phi * dir * magCellVel);
	}
    }
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::sbPartialSlipMooneyFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "shaftOmega", omega_);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "m", m_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sbPartialSlipMooneyFvPatchVectorField
    );
}
