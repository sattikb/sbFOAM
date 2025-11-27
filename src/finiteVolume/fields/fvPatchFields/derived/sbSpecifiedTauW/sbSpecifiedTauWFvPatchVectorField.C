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
    fixedGradientFvPatchVectorField(p, iF),
    f_()
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    f_(dict.lookup<scalar>("f"))
{
    // Initialize with zero gradient, will be updated in updateCoeffs
    gradient() = vectorField(p.size(), vector::zero);
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
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
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
    f_(ptf.f_)
{}

sbSpecifiedTauWFvPatchVectorField::sbSpecifiedTauWFvPatchVectorField
(
    const sbSpecifiedTauWFvPatchVectorField& sbstw,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(sbstw, iF),
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


    vectorField& dudn = gradient();
    forAll(dudn, faceI)
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
	scalar tauLim = stressWN * f_;
	
	dudn[faceI] = gradUf & nf;
	if (tauW > tauLim) 
	{
		//dudn[faceI] = f_ * pp[faceI] * vector::zero;
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

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void sbSpecifiedTauWFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    writeEntry(os, "f", f_);
}

makePatchTypeField
(
    fvPatchVectorField,
    sbSpecifiedTauWFvPatchVectorField
);

} // End namespace Foam
