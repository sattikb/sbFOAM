#include "rotatingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatingWallVelocityFvPatchVectorField::
rotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(),
    axis_(Zero),
    omega_(0)
{}


Foam::rotatingWallVelocityFvPatchVectorField::
rotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(Function1<scalar>::New("omega", dict))
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


Foam::rotatingWallVelocityFvPatchVectorField::
rotatingWallVelocityFvPatchVectorField
(
    const rotatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_, false)
{}


Foam::rotatingWallVelocityFvPatchVectorField::
rotatingWallVelocityFvPatchVectorField
(
    const rotatingWallVelocityFvPatchVectorField& rwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(rwvpvf, iF),
    origin_(rwvpvf.origin_),
    axis_(rwvpvf.axis_),
    omega_(rwvpvf.omega_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    scalar om = omega_->value(t);

    // Calculate the rotating wall velocity from the specification of the motion
    const vectorField Up
    (
        (-om)*((patch().Cf() - origin_) ^ (axis_/mag(axis_)))
    );

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    const vectorField n(patch().nf());
    vectorField::operator=(Up - n*(n & Up));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::rotatingWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "axis", axis_);
    writeEntry(os, omega_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rotatingWallVelocityFvPatchVectorField
    );
}
