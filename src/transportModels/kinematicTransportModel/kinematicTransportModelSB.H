// BASE CLASS FOR ALL TRANSPORT MODELS USED BY INCOMPRESSIBLE TURBULENCE 
// MODELS

#ifndef kinematicTransportModelSB_H
#define kinematicTransportModelSB_H

#include "primitiveFieldsFwdSB.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class kinematicTransportModel
{
public:

    //- Runtime type information
    TypeName("kinematicTransportModel");


    // Constructors

        //- Construct from components
        kinematicTransportModel();

        //- Disallow default bitwise copy construction
        kinematicTransportModel(const kinematicTransportModel&);


    //- Destructor
    virtual ~kinematicTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const = 0;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kinematicTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
