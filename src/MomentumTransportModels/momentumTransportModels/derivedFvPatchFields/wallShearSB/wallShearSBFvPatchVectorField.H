// Description
//     Set a constant shear stress as tau0 = -nuEff dU/dn.

#ifndef wallShearSBFvPatchVectorField_H
#define wallShearSBFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class wallShearSBFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class wallShearSBFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Constant shear stress
        const vector tau0_;


public:

    //- Runtime type information
    TypeName("wallShearSB");


    // Constructors

        //- Construct from patch and internal field
        wallShearSBFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallShearSBFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        wallShearSBFvPatchVectorField
        (
            const wallShearSBFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        wallShearSBFvPatchVectorField
        (
            const wallShearSBFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        wallShearSBFvPatchVectorField
        (
            const wallShearSBFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new wallShearSBFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
