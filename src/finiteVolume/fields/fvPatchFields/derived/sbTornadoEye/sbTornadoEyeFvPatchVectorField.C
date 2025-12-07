#include "sbTornadoEyeFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sbTornadoEyeFvPatchVectorField::sbTornadoEyeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    omega_(),
    eyeLoc_(),
    rLim_()
{}


Foam::sbTornadoEyeFvPatchVectorField::sbTornadoEyeFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    omega_(dict.lookup<vector>("shaftOmega")),
    eyeLoc_(dict.lookup<vector>("eyeLoc")),
    rLim_(dict.lookup<scalar>("rLim"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::sbTornadoEyeFvPatchVectorField::sbTornadoEyeFvPatchVectorField
(
    const sbTornadoEyeFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    omega_(ptf.omega_),
    eyeLoc_(ptf.eyeLoc_),
    rLim_(ptf.rLim_)
{}


Foam::sbTornadoEyeFvPatchVectorField::sbTornadoEyeFvPatchVectorField
(
    const sbTornadoEyeFvPatchVectorField& pstauSB,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pstauSB, iF),
    omega_(pstauSB.omega_),
    eyeLoc_(pstauSB.eyeLoc_),
    rLim_(pstauSB.rLim_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sbTornadoEyeFvPatchVectorField::updateCoeffs()
{
    Info<<"SATTIK IN MOONEY PARTIAL SLIP."<<endl;
    if (updated())
    {
        return;
    }

    const fvPatch& p = patch();
    const vectorField n = p.nf();

    vectorField& Up = *this;
    forAll(Up, faceI)
    {
	const vector rad = patch().Cf()[faceI] - eyeLoc_;
	
	Up[faceI] = vector::zero;
	
	if (mag(rad)<=rLim_)
	{
		const vector plateVel = omega_ ^ rad;
		Up[faceI] = plateVel;
	}
    }
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::sbTornadoEyeFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "shaftOmega", omega_);
    writeEntry(os, "eyeLoc", eyeLoc_);
    writeEntry(os, "rLim", rLim_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        sbTornadoEyeFvPatchVectorField
    );
}
