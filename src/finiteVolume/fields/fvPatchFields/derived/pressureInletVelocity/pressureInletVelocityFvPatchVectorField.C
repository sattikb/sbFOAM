#include "pressureInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
    const pressureInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
    const pressureInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    tmp<vectorField> n = patch().nf();
    const Field<scalar>& magS = patch().magSf();

    if (phi.dimensions() == dimFlux)
    {
        operator==(n*phip/magS);
    }
    else if (phi.dimensions() == dimMassFlux)
    {
	    Info<<"SATTIK P-U BC"<<endl;
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(n*phip/(rhop*magS));
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pressureInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressureInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
