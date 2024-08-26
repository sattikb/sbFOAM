// BASE CLASS FOR DERIVED THERMO PROPERTIES BASED ON DENSITY

#ifndef sbRhoThermo_H
#define sbRhoThermo_H

#include "fluidThermoSB.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class sbRhoThermo Declaration
\*---------------------------------------------------------------------------*/

class sbRhoThermo
:
    virtual public fluidThermo
{
public:

    // Public Classes

        //- Forward declare the implementation class
        class implementation;

        //- Forward declare the composite class
        class composite;


    //- Runtime type information
    TypeName("sbRhoThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        sbRhoThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<sbRhoThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~sbRhoThermo();


    // Member Functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho() = 0;

            //- Old-time density [kg/m^3]
            virtual tmp<volScalarField> rho0() const = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho) = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;
};


/*---------------------------------------------------------------------------*\
                  Class sbRhoThermo::implementation Declaration
\*---------------------------------------------------------------------------*/

class sbRhoThermo::implementation
:
    virtual public sbRhoThermo
{
protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'thermo:rho' to avoid (potential) conflict with solver density
        volScalarField rho_;

        //- Compressibility [s^2/m^2]
        volScalarField psi_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;


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

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Old-time density [kg/m^3]
            virtual tmp<volScalarField> rho0() const;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho);

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implementation&) = delete;
};


/*---------------------------------------------------------------------------*\
                    Class sbRhoThermo::composite Declaration
\*---------------------------------------------------------------------------*/

class sbRhoThermo::composite
:
    public basicThermo::implementation,
    public fluidThermo::implementation,
    public sbRhoThermo::implementation
{
public:

    // Constructors

        //- Construct from mesh and phase name
        composite
        (
            const fvMesh& mesh,
            const word& phaseName
        )
        :
            basicThermo::implementation(mesh, phaseName),
            fluidThermo::implementation(mesh, phaseName),
            sbRhoThermo::implementation(mesh, phaseName)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
