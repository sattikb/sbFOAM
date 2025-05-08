#ifndef tauRTFvPatchVectorField_H
#define tauRTFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tauRTFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data
	scalar tauRT_;
	scalar mu_;

public:

    //- Runtime type information
    TypeName("tauRT_SB");


    // Constructors

        //- Construct from patch and internal field
        tauRTFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tauRTFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rotatingWallVelocityFvPatchVectorField
        //  onto a new patch
        tauRTFvPatchVectorField
        (
            const tauRTFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        tauRTFvPatchVectorField
        (
            const tauRTFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        tauRTFvPatchVectorField
        (
            const tauRTFvPatchVectorField&,
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
                new tauRTFvPatchVectorField(*this, iF)
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

#endif
