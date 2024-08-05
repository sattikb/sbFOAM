// THIS IS THE BASE CLASS FOR ALL THE TRANSPORT MODELS USED BY THE 
// COMPRESSIBLE TURBULENCE MODELS.

#ifndef dynamicTransportModelSB_H
#define dynamicTransportModelSB_H

#include "scalarField.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dynamicTransportModel Declaration
\*---------------------------------------------------------------------------*/

class dynamicTransportModel
{
public:

    //- Runtime type information
    TypeName("dynamicTransportModel");


    // Constructors

        //- Construct from components
        dynamicTransportModel();

        //- Disallow default bitwise copy construction
        dynamicTransportModel
        (
            const dynamicTransportModel&
        ) = delete;


    //- Destructor
    virtual ~dynamicTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> mu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
