// Description
//     Set a constant shear stress as tau0 = -nuEff dU/dn.

#ifndef fixedShearStressFvPatchVectorField_H
#define fixedShearStressFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class fixedShearStressFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedShearStressFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Constant shear stress
        const vector tau0_;


public:

    //- Runtime type information
    TypeName("fixedShearStress");


    // Constructors

        //- Construct from patch and internal field
        fixedShearStressFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedShearStressFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        fixedShearStressFvPatchVectorField
        (
            const fixedShearStressFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fixedShearStressFvPatchVectorField
        (
            const fixedShearStressFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        fixedShearStressFvPatchVectorField
        (
            const fixedShearStressFvPatchVectorField&,
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
                new fixedShearStressFvPatchVectorField(*this, iF)
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
