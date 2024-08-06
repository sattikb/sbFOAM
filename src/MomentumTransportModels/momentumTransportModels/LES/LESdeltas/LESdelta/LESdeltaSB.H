// ABSTRACT BASE CLASS FOR LES DELTAS

#ifndef LESdeltaSB_H
#define LESdeltaSB_H

#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class LESdelta
{
protected:

    // Protected data

        const momentumTransportModel& momentumTransportModel_;

        volScalarField delta_;


public:

    //- Runtime type information
    TypeName("LESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESdelta,
            dictionary,
            (
                const word& name,
                const momentumTransportModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, momentumTransportModel and dictionary
        LESdelta
        (
            const word& name,
            const momentumTransportModel& turbulence
        );

        //- Disallow default bitwise copy construction
        LESdelta(const LESdelta&) = delete;


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary& dict
        );

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const momentumTransportModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable&
        );


    //- Destructor
    virtual ~LESdelta()
    {}


    // Member Functions

        //- Return momentumTransportModel reference
        const momentumTransportModel& turbulence() const
        {
            return momentumTransportModel_;
        }

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        void operator=(const LESdelta&) = delete;

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
