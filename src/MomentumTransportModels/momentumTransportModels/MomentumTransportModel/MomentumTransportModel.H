// TEMPLATED ABSTRACT BASE CLASS FOR TURBULENCE MODELS

#ifndef MomentumTransportModel_H
#define MomentumTransportModel_H

#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Alpha,
    class Rho,
    class BasicMomentumTransportModel,
    class TransportModel
>
class MomentumTransportModel
:
    public BasicMomentumTransportModel
{

public:

    typedef Alpha alphaField;
    typedef Rho rhoField;
    typedef TransportModel transportModel;


protected:

    // Protected data

        const alphaField& alpha_;
        const transportModel& transport_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            MomentumTransportModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const transportModel& transport
            ),
            (alpha, rho, U, alphaRhoPhi, phi, transport)
        );


    // Constructors

        //- Construct
        MomentumTransportModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );

        //- Disallow default bitwise copy construction
        MomentumTransportModel(const MomentumTransportModel&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<MomentumTransportModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    //- Destructor
    virtual ~MomentumTransportModel()
    {}


    // Member Functions

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Access function to incompressible transport model
        const transportModel& transport() const
        {
            return transport_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MomentumTransportModel.C"
#endif

#endif
