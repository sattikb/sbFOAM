//TEMPLATE ABSTRACT BASE CLASS FOR SINGLE-PHASE COMPRESSIBLE
//TURBULENCE MODELS

#ifndef CompressibleMomentumTransportModelSB_H
#define CompressibleMomentumTransportModelSB_H

#include "MomentumTransportModelSB.H"
#include "compressibleMomentumTransportModelSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class CompressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class TransportModel>
class CompressibleMomentumTransportModel
:
    public MomentumTransportModel
    <
        geometricOneField,
        volScalarField,
        compressibleMomentumTransportModel,
        TransportModel
    >
{

public:

    typedef geometricOneField alphaField;
    typedef volScalarField rhoField;
    typedef TransportModel transportModel;


    // Constructors

        //- Construct
        CompressibleMomentumTransportModel
        (
            const word& type,
            const geometricOneField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<CompressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& transportModel
        );


    //- Destructor
    virtual ~CompressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const
        {
            return this->transport_.mu();
        }

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return this->transport_.mu(patchi);
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return this->transport_.mu()/this->rho_;
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return
                this->transport_.mu(patchi)
               /this->rho_.boundaryField()[patchi];
        }

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const
        {
            return this->rho_*this->nut();
        }

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const
        {
            return this->rho_.boundaryField()[patchi]*this->nut(patchi);
        }

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return mut() + mu();
        }

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const
        {
            return mut(patchi) + mu(patchi);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompressibleMomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
