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
    f_(dict.lookup<scalar>("f"))
{
    this->refValue() = vector::zero;
    this->refGrad()  = vector::zero;
    this->valueFraction() = 1.0;
    // Initialize with zero gradient, will be updated in updateCoeffs
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
    f_(ptf.f_)
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const sbSpecifiedTauWFvPatchVectorField& sbstw,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(sbstw, iF),
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

//    const fvPatchVectorField& Ubound = U.boundaryField()[p.index()];
    const fvPatchScalarField& pp = pField.boundaryField()[p.index()];
    const fvPatchScalarField& mup = mu.boundaryField()[p.index()];

    tensorField gradU = fvc::grad(U)().boundaryField()[p.index()];

    // Compute gradient: du/dn = f * p


//    vectorField& Up = *this;
    this->refValue() = vector::zero;
    this->refGrad()  = vector::zero;
    this->valueFraction() = 1.0;

    forAll(n, faceI)
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
	scalar tauLim = -stressWN * f_;
	
	if (tauW < tauLim) 
	{
		this->valueFraction()[faceI] = 1.0;
	}
	else
	{
		this->valueFraction()[faceI] = 0.0;
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
    writeEntry(os, "f", f_);
}

makePatchTypeField
(
    fvPatchVectorField,
    sbSpecifiedTauWFvPatchVectorField
);

} // End namespace Foam
