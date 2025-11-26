#include "sbPartialSlipOnlyTauFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sbPartialSlipOnlyTauFvPatchVectorField::sbPartialSlipOnlyTauFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    movingWallVelocityFvPatchVectorField(p, iF),// fixedValueFvPatchVectorField(p, iF),
    omega_(),
    alpha_(),
    m_()
{}


Foam::sbPartialSlipOnlyTauFvPatchVectorField::sbPartialSlipOnlyTauFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    movingWallVelocityFvPatchVectorField(p, iF), //fixedValueFvPatchVectorField(p, iF, dict, false),
    omega_(dict.lookup<vector>("shaftOmega")),
    alpha_(dict.lookup<scalar>("alpha")),
    m_(dict.lookup<scalar>("m"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::sbPartialSlipOnlyTauFvPatchVectorField::sbPartialSlipOnlyTauFvPatchVectorField
(
    const sbPartialSlipOnlyTauFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    movingWallVelocityFvPatchVectorField(ptf, p, iF, mapper), // fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_),
    alpha_(ptf.alpha_),
    m_(ptf.m_)
{}


Foam::sbPartialSlipOnlyTauFvPatchVectorField::sbPartialSlipOnlyTauFvPatchVectorField
(
    const sbPartialSlipOnlyTauFvPatchVectorField& pstauSB,
    const DimensionedField<vector, volMesh>& iF
)
:
    movingWallVelocityFvPatchVectorField(pstauSB, iF),// fixedValueFvPatchVectorField(pstauSB, iF),
    omega_(pstauSB.omega_),
    alpha_(pstauSB.alpha_),
    m_(pstauSB.m_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sbPartialSlipOnlyTauFvPatchVectorField::updateCoeffs()
{
    Info<<"SATTIK IN NORMALIZED MOONEY PARTIAL SLIP."<<endl;

    Info << "----- patch: " << patch().name()
	 << "  type: " << patch().type()
	 << "  local faces: " << patch().size()
	 << "  global faces: " << returnReduce(patch().size(), sumOp<label>())
	 << endl;
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

//    const fvPatchVectorField& Ubound = U.boundaryField()[p.index()];
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
	scalar tauStar = - tauW / stressWN;

	const vector plateVel = omega_ ^ patch().Cf()[faceI];

	Info<<"------ value of taustar is: "<<tauStar<<endl;

        if (tauStar>0.2)
        {
		localSlipFaces++;
	    vector dir = -stressWS/tauW;
            scalar coeff = pow(alpha_*tauStar, m_);

	    label own = p.faceCells()[faceI];
            vector cellVel = U[own];
            scalar magCellVel = mag(cellVel);

	Info<<"------ value of vs is: "<<coeff<<" and chi is: "<<tauStar<<endl;
            Up[faceI] = plateVel + (coeff * dir * magCellVel);
        }
        else
        {
            Up[faceI] = plateVel;
        }
    }
    
    label globalProcessed = returnReduce(localProcessedFaces, sumOp<label>());
    label globalSlip      = returnReduce(localSlipFaces,      sumOp<label>());

    Info << "=== MOONEY SLIP BC ACTIVE === "
         << "processed " << globalProcessed
         << " faces globally, slip on " << globalSlip
	 << " faces this update (patch: " << patch().name() << ")"
	 << endl;
    //fixedValueFvPatchVectorField::updateCoeffs();
    movingWallVelocityFvPatchVectorField::updateCoeffs();
}


void Foam::sbPartialSlipOnlyTauFvPatchVectorField::write(Ostream& os) const
{
    movingWallVelocityFvPatchVectorField::write(os);  // fvPatchVectorField::write(os);
    writeEntry(os, "shaftOmega", omega_);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "m", m_);
    writeEntry(os, "value", *this);
}

bool Foam::sbPartialSlipOnlyTauFvPatchVectorField::moving() const
{
	    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sbPartialSlipOnlyTauFvPatchVectorField
    );
}
