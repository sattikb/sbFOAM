#ifndef sbTornadoEyeFvPatchVectorField_H
#define sbTornadoEyeFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class sbTornadoEyeFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sbTornadoEyeFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    // Private Data
    vector omega_; 
    vector eyeLoc_;
    scalar rLim_;    

public:
    //- Runtime type information
    TypeName("sbTornadoEye");

    // Constructors
    //- Construct from patch and internal field
    sbTornadoEyeFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    //- Construct from patch, internal field, and dictionary
    sbTornadoEyeFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given sbTornadoEyeFvPatchVectorField onto a new patch
    sbTornadoEyeFvPatchVectorField
    (
        const sbTornadoEyeFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Copy constructor
    sbTornadoEyeFvPatchVectorField
    (
        const sbTornadoEyeFvPatchVectorField&
    ) = delete;

    //- Copy constructor setting internal field reference
    sbTornadoEyeFvPatchVectorField
    (
        const sbTornadoEyeFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    //- Create a clone with internal field reference
    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new sbTornadoEyeFvPatchVectorField(*this, iF)
        );
    }

    // Member Functions
    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif
