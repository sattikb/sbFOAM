#ifndef tauRTWallFvPatchVectorField_H
#define tauRTWallFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rotatingWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class tauRTWallFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data
	scalar tauRT_;
	scalar rIn_;
	scalar rOut_;
	scalar omegaIn_;

public:

    //- Runtime type information
    TypeName("tauRT_FullDeveloped_SB");


    // Constructors

        //- Construct from patch and internal field
        tauRTWallFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tauRTWallFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rotatingWallVelocityFvPatchVectorField
        //  onto a new patch
        tauRTWallFvPatchVectorField
        (
            const tauRTWallFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        tauRTWallFvPatchVectorField
        (
            const tauRTWallFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        tauRTWallFvPatchVectorField
        (
            const tauRTWallFvPatchVectorField&,
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
                new tauRTWallFvPatchVectorField(*this, iF)
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
