#include "tauRTWallFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::tauRTWallFvPatchVectorField::
tauRTWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    tauRT_(),
    rIn_(),
    rOut_(),
    omegaIn_()
{}


Foam::tauRTWallFvPatchVectorField::
tauRTWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    tauRT_(dict.lookupOrDefault<scalar>("tauRT", 0.0)),
    rIn_(dict.lookupOrDefault<scalar>("rIn", 0.0)),
    rOut_(dict.lookupOrDefault<scalar>("rOut", 0.0)),
    omegaIn_(dict.lookupOrDefault<scalar>("omegaIn",0.0))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // Evaluate the wall velocity
        updateCoeffs();
    }
}


Foam::tauRTWallFvPatchVectorField::
tauRTWallFvPatchVectorField
(
    const tauRTWallFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    tauRT_(ptf.tauRT_),
    rIn_(ptf.rIn_),
    rOut_(ptf.rOut_),
    omegaIn_(ptf.omegaIn_)
{}


Foam::tauRTWallFvPatchVectorField::
tauRTWallFvPatchVectorField
(
    const tauRTWallFvPatchVectorField& trtFDpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(trtFDpvf, iF),
    tauRT_(trtFDpvf.tauRT_),
    rIn_(trtFDpvf.rIn_),
    rOut_(trtFDpvf.rOut_),
    omegaIn_(trtFDpvf.omegaIn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tauRTWallFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//    const volScalarField& muSBField = db().lookupObject<volScalarField>("thermo:mu");
//    scalar muSB = average(muSBField.internalField()).value();
    const scalar num1 = rOut_*rOut_ - rIn_*rIn_;
    const scalar den1 = rIn_*rIn_ * rOut_*rOut_;
    const scalar num2 = rOut_*rOut_ * tauRT_;
    const scalar den2 = 2*scalar(50000);
    const scalar term1 = num1/den1;
    const scalar term2 = num2/den2;
    const scalar omegaOut = omegaIn_ + term1*term2;

    // Calculate the rotating wall velocity from the specification of the motion
    const vectorField Up
    (
       // (-omegaOut)*( patch().Cf() ^ vector(0,0,1) )
       vector(0,0,omegaOut) ^ patch().Cf()
    );

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    const vectorField n(patch().nf());
    vectorField::operator=(Up - n*(n & Up));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::tauRTWallFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "tauRT", tauRT_);
    writeEntry(os, "rIn", rIn_);
    writeEntry(os, "rOut", rOut_);
    writeEntry(os, "omegaIn", omegaIn_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tauRTWallFvPatchVectorField
    );
}
