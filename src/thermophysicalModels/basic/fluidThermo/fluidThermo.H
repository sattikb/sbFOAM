// ABC FOR FLUIDTHERMO, WHICH IS USED IN THE FLUIDTHERMO::IMPLEMENTATION
// CLASS
// Foam::fluidThermo

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo.H"
#include "kinematicTransportModelSB.H"	//FOR INCOMPRESSIBLE TRANSPORT MODELS
#include "dynamicTransportModelSB.H"	//FOR COMPRESSIBLE TRANSPORT MODELS

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    virtual public basicThermo,
    public kinematicTransportModel,
    public dynamicTransportModel
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;


    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<fluidThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~fluidThermo();


    // Member Functions

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p() = 0;

            //- Pressure [Pa]
            virtual const volScalarField& p() const = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho) = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Fields derived from thermodynamic state variables

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const = 0;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;


        // Fields derived from transport state variables

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const = 0;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;
};


/*---------------------------------------------------------------------------*\
                 Class fluidThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo::implementation
:
    virtual public fluidThermo
{
protected:

    // Protected data

        // Fields

            //- Pressure [Pa]
            volScalarField& p_;


public:

    // Constructors

        //- Construct from mesh and phase name
        implementation(const fvMesh&, const word& phaseName);

        //- Construct from mesh, dictionary and phase name
        implementation(const fvMesh&, const dictionary&, const word& phaseName);

        //- Disallow default bitwise copy construction
        implementation(const implementation&) = delete;


    //- Destructor
    virtual ~implementation();


    // Member Functions

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
