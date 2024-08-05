// TEMPLATED ABSTRACT BASE CLASS FOR TURBULENCE MODELS

#ifndef compressibleMomentumTransportModelSB_H
#define compressibleMomentumTransportModelSB_H

#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
             Class compressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleMomentumTransportModel
:
    public momentumTransportModel
{

protected:

    // Protected data

        const volScalarField& rho_;


public:

    //- Runtime type information
    TypeName(momentumTransportModel::typeName_());


    // Constructors

        //- Construct from components
        compressibleMomentumTransportModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        compressibleMomentumTransportModel
        (
            const compressibleMomentumTransportModel&
        ) = delete;


    //- Destructor
    virtual ~compressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Return the density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devTau() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleMomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
