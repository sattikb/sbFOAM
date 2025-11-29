#include "sbSpecifiedTauWFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"

namespace Foam
{

// Constructors
sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    omega_(),
    f_()
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    omega_(dict.lookup<vector>("shaftOmega")),
    f_(dict.lookup<scalar>("f"))
{
    this->refValue() = vector::zero;
    this->refGrad()  = vector::zero;
    this->valueFraction() = 1.0;

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
	this->refValue() = *this;
        this->valueFraction() = 1.0;
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }
}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const sbSpecifiedTauWFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_),
    f_(ptf.f_)
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const sbSpecifiedTauWFvPatchVectorField& sbstw,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(sbstw, iF),
    omega_(sbstw.omega_),
    f_(sbstw.f_)
{}

// Member Functions
void sbSpecifiedTauWFvPatchVectorField::updateCoeffs()
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
    // Get viscosity, pressure and velocity field
    const volScalarField& pField = db().lookupObject<volScalarField>("p");
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    const volScalarField& mu = db().lookupObject<volScalarField>("thermo:mu");

    const fvPatchScalarField& pp = pField.boundaryField()[p.index()];
    const fvPatchScalarField& mup = mu.boundaryField()[p.index()];

    tensorField gradU = fvc::grad(U)().boundaryField()[p.index()];

//    u_face = (1-phi)(u_cc + dUdnp*ds) + phi*Up
    vectorField& Up = this->refValue();
    vectorField& dUdnp = this->refGrad();
    scalarField& phi = this->valueFraction();

    Up = omega_ ^ p.Cf();
    dUdnp = vector::zero;
    phi = 1.0;

    forAll(n, faceI)
    {
	localProcessedFaces++;

        scalar muf = mup[faceI];
        scalar pf = pp[faceI];
	const vector plateVel = omega_ ^ p.Cf()[faceI];

        const vector& nf = n[faceI];
        const tensor& gradUf = gradU[faceI];

	tensor tauf = muf*(gradUf + gradUf.T()) - (2.0/3.0)*muf*(tr(gradUf))*I;
	tensor sigmaf = tauf - pf*I;

	vector traction = sigmaf & nf;
	scalar stressWN = traction & nf;
	vector stressWS = traction  - stressWN*nf;
	
	scalar tauW = mag(stressWS);
	scalar tauLim = -stressWN * f_;
	
	Info<<"The wall shear stress is: "<<tauW<<" and the tauLim is : "<<tauLim<<endl;
	if (tauW < tauLim) 
	{
		phi[faceI] = 1.0;
		Up[faceI] = plateVel;
	}
	else
	{
		phi[faceI] = 0.0;
		dUdnp[faceI] = vector::zero;
		localSlipFaces++;
	}
    }

    label globalProcessed = returnReduce(localProcessedFaces, sumOp<label>());
    label globalSlip      = returnReduce(localSlipFaces,      sumOp<label>());
    Info << "=== COULOMB SLIP BC ACTIVE === "
         << "processed " << globalProcessed
         << " faces globally, slip on " << globalSlip
	 << " faces this update (patch: " << patch().name() << ")"
	 << endl;

    mixedFvPatchVectorField::updateCoeffs();
}

void sbSpecifiedTauWFvPatchVectorField::write(Ostream& os) const
{
    mixedFvPatchVectorField::write(os);
    writeEntry(os, "shaftOmega", omega_);
    writeEntry(os, "f", f_);
}

makePatchTypeField
(
    fvPatchVectorField,
    sbSpecifiedTauWFvPatchVectorField
);

} // End namespace Foam
