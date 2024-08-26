// WNERGY BASED ON DENSITY, BUT WITH SOME CASE MODIFICATIONS

#ifndef sbThermo_H
#define sbThermo_H

#include "sbRhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sbThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicRhoThermo, class MixtureType>
class sbThermo
:
    public heThermo<BasicRhoThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();


public:

    //- Runtime type information
    TypeName("sbThermo");


    // Constructors

        //- Construct from mesh and phase name
        sbThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        sbThermo(const sbThermo<BasicRhoThermo, MixtureType>&) = delete;


    //- Destructor
    virtual ~sbThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const sbThermo<BasicRhoThermo, MixtureType>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "sbThermo.C"
#endif

#endif
