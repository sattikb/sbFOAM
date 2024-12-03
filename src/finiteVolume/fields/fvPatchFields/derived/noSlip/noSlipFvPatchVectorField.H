//    This boundary condition fixes the velocity to zero at walls.
//
//    Example of the boundary condition specification:
//    \verbatim
//    <patchName>
//    {
//        type            noSlip;
//    }
//    \endverbatim

#ifndef noSlipFvPatchVectorField_H
#define noSlipFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class noSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class noSlipFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("noSlip");


    // Constructors

        //- Construct from patch and internal field
        noSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        noSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given noSlipFvPatchVectorField
        //  onto a new patch
        noSlipFvPatchVectorField
        (
            const noSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        noSlipFvPatchVectorField
        (
            const noSlipFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        noSlipFvPatchVectorField
        (
            const noSlipFvPatchVectorField&,
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
                new noSlipFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
