#include "tauRTFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::tauRTFvPatchVectorField::
tauRTFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    tauRT_(),
    mu_()
{}


Foam::tauRTFvPatchVectorField::
tauRTFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    tauRT_(dict.lookupOrDefault<scalar>("tauRT", 0.0)),
    mu_(dict.lookupOrDefault<scalar>("mu",0.0))
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


Foam::tauRTFvPatchVectorField::
tauRTFvPatchVectorField
(
    const tauRTFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    tauRT_(ptf.tauRT_),
    mu_(ptf.mu_)
{}


Foam::tauRTFvPatchVectorField::
tauRTFvPatchVectorField
(
    const tauRTFvPatchVectorField& trtpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(trtpvf, iF),
    tauRT_(trtpvf.tauRT_),
    mu_(trtpvf.mu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tauRTFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

        const fvMesh& mesh = this->patch().boundaryMesh().mesh();
        const labelList& faceCells = this->patch().faceCells();
        const vectorField& faceCenters = this->patch().Cf();

	const volVectorField& U = db().lookupObject<volVectorField>(this->internalField().name());
        Field<vector>& patchField = *this;

        forAll(patchField, facei)
        {
            // Get internal cell data
            const label celli = faceCells[facei];

	    const vector& U_cell = U.internalField()[celli];
            
	    const vector& faceCenter = faceCenters[facei];
	    const vector& cellCenter = mesh.C()[celli];

	    const vector d = faceCenter - cellCenter;
	    scalar dx = d[0];
	    scalar dy = d[1];

            // Compute cylindrical coordinates
            const scalar xf = faceCenter.x();
            const scalar yf = faceCenter.y();
            const scalar rf = sqrt( xf*xf + yf*yf );

	    scalar a1 = xf/rf;
	    scalar b1 = yf/rf;
	    scalar c1 = 0.0;

	    
	    scalar a2a =   yf;
	    scalar a2b = - xf * yf / dx;
	    scalar a2c = - yf * yf / dy;
	    scalar b2a =   xf;
	    scalar b2b = - xf * xf / dx;
	    scalar b2c = - xf * yf / dy;
	    scalar c2a =   tauRT_ / mu_ * rf*rf;
	    scalar c2b = - ( yf*yf/dy + xf*yf/dx ) * U_cell[0];
	    scalar c2c =   ( xf*xf/dx + xf*yf/dy ) * U_cell[1];

	    scalar a2 =   ( a2a + a2b + a2c );
	    scalar b2 = - ( b2a + b2b + b2c );
	    scalar c2 =   ( c2a + c2b + c2c );

	    scalar det = a1*b2 - a2*b1;
	    
	    //Cramer's Rule
	    scalar ux = ( c1*b2 - c2*b1) / det;
	    scalar uy = (-c1*a2 + c2*a1) / det;
	    scalar uz = 0.0;

            vector U_p(ux, uy, uz);

            patchField[facei] = U_p;
        }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::tauRTFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "tauRT", tauRT_);
    writeEntry(os, "mu", mu_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tauRTFvPatchVectorField
    );
}
