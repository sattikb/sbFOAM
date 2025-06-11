#include "shearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "kinematicMomentumTransportModel.H"
#include "dynamicMomentumTransportModelSB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(shearStress, 0);
    addToRunTimeSelectionTable(functionObject, shearStress, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::shearStress::shearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, false),
    phaseName_(word::null)
{
    read(dict);
    resetLocalObjectName(IOobject::groupName(type(), phaseName_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::shearStress::~shearStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::shearStress::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


bool Foam::functionObjects::shearStress::execute()
{
    const word fieldName(IOobject::groupName(type(), phaseName_));

    typedef compressibleMomentumTransportModel cmpModel;
    typedef incompressibleMomentumTransportModel icoModel;

    const word momentumTransportModelName
    (
        IOobject::groupName(momentumTransportModel::typeName, phaseName_)
    );

    if (mesh_.foundObject<cmpModel>(momentumTransportModelName))
    {
        const cmpModel& model =
            mesh_.lookupObject<cmpModel>(momentumTransportModelName);

	Info<<"SATTIK returned devTau"<<endl;
        return store(fieldName, model.devTau());
    }
    else if (mesh_.foundObject<icoModel>(momentumTransportModelName))
    {
        const icoModel& model =
            mesh_.lookupObject<icoModel>(momentumTransportModelName);

	Info<<"SATTIK returned devSigma"<<endl;
        return store(fieldName, model.devSigma());
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model "
            << momentumTransportModelName << " in the database"
            << exit(FatalError);

        return false;
    }
}


bool Foam::functionObjects::shearStress::write()
{
    return writeLocalObjects::write();
}
