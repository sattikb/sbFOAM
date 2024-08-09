// ABC FOR THERMOPHYSICAL TRANSPORT MODELS (LAMINAR, RAS, LES)

#ifndef thermophysicalTransportModelSB_H
#define thermophysicalTransportModelSB_H

#include "compressibleMomentumTransportModelSB.H"
#include "fluidThermoSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class thermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class thermophysicalTransportModel
:
    public IOdictionary
{
protected:

    // Protected data

        const compressibleMomentumTransportModel& momentumTransportModel_;


public:

    //- Runtime type information
    TypeName("thermophysicalTransport");


    // Constructors

        //- Construct from compressibleMomentumTransportModel
        thermophysicalTransportModel
        (
            const compressibleMomentumTransportModel& momentumTransport
        );

        //- Disallow default bitwise copy construction
        thermophysicalTransportModel
        (
            const thermophysicalTransportModel&
        ) = delete;


    //- Destructor
    virtual ~thermophysicalTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        const compressibleMomentumTransportModel& momentumTransport() const
        {
            return momentumTransportModel_;
        }

        //- Access function to incompressible transport model
        virtual const fluidThermo& thermo() const = 0;

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const = 0;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const = 0;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const = 0;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const = 0;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const = 0;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const = 0;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const = 0;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const = 0;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const = 0;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const = 0;

        //- Solve the thermophysical transport model equations
        //  and correct the transport coefficients
        virtual void correct() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
