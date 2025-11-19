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
    fixedValueFvPatchVectorField(p, iF),
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
    fixedValueFvPatchVectorField(p, iF, dict, false),
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
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
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
    fixedValueFvPatchVectorField(pstauSB, iF),
    omega_(pstauSB.omega_),
    alpha_(pstauSB.alpha_),
    m_(pstauSB.m_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sbPartialSlipOnlyTauFvPatchVectorField::updateCoeffs()
{
    Info<<"SATTIK IN MODIFIED MOONEY PARTIAL SLIP."<<endl;
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


        if (tauStar>0.2)
        {
	    vector dir = -stressWS/tauW;
            scalar coeff = pow(alpha_*tauStar, m_);

	    label own = p.faceCells()[faceI];
            vector cellVel = U[own];
            scalar magCellVel = mag(cellVel);

            Up[faceI] = plateVel + (coeff * dir * magCellVel);
        }
        else
        {
            Up[faceI] = plateVel;
        }
//	Info<<"------ value of taustar is: "<<tauStar<<endl;
    }
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::sbPartialSlipOnlyTauFvPatchVectorField::write(Ostream& os) const
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
        sbPartialSlipOnlyTauFvPatchVectorField
    );
}
