#include "thermophysicalTransportModelSB.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermophysicalTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermophysicalTransportModel::thermophysicalTransportModel
(
    const compressibleMomentumTransportModel& momentumTransport
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName
            (
                typeName, momentumTransport.alphaRhoPhi().group()
            ),
            momentumTransport.time().constant(),
            momentumTransport.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),

    momentumTransportModel_(momentumTransport)
{
    // Add run-time re-reading of thermophysicalTransport dictionary
    // after construction to avoid problems if the dictionary is not present
    readOpt() = IOobject::MUST_READ_IF_MODIFIED;
    addWatch();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::thermophysicalTransportModel::read()
{
    return regIOobject::read();
}


void Foam::thermophysicalTransportModel::correct()
{}


// ************************************************************************* //
