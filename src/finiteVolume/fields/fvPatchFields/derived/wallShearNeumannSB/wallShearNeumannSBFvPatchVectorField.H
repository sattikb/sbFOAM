#ifndef wallShearNeumannSBFvPatchVectorField_H
#define wallShearNeumannSBFvPatchVectorField_H

#include "fixedGradientFvPatchFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class wallShearNeumannSBFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private Data

    vectorField traction_;  // Traction vector (input)
    scalar mu_;             // Dynamic viscosity (input)


public:

    //- Runtime type information
    TypeName("wallShearNeumannSB");


    // Constructors

        //- Construct from patch and internal field
        wallShearNeumannSBFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field, and dictionary
        wallShearNeumannSBFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping
        wallShearNeumannSBFvPatchVectorField
        (
            const wallShearNeumannSBFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wallShearNeumannSBFvPatchVectorField
        (
            const wallShearNeumannSBFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new wallShearNeumannSBFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        wallShearNeumannSBFvPatchVectorField
        (
            const wallShearNeumannSBFvPatchVectorField&,
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
                new wallShearNeumannSBFvPatchVectorField(*this, iF)
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
