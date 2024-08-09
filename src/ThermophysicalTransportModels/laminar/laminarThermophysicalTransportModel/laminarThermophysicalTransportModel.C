#include "laminarThermophysicalTransportModelSB.H"
#include "unityLewisFourier.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
void Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::printCoeffs
(
    const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::laminarThermophysicalTransportModel
(
    const word& type,
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    BasicThermophysicalTransportModel(momentumTransport, thermo),
    laminarDict_(this->subOrEmptyDict("laminar")),
    printCoeffs_(laminarDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(laminarDict_.optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
Foam::autoPtr
<
    Foam::laminarThermophysicalTransportModel
    <
        BasicThermophysicalTransportModel
    >
>
Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::New
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
{
    IOobject header
    (
        IOobject::groupName
        (
            thermophysicalTransportModel::typeName,
            momentumTransport.alphaRhoPhi().group()
        ),
        momentumTransport.time().constant(),
        momentumTransport.mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (header.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary modelDict(header);

        const word modelType(modelDict.subDict("laminar").lookup( "model"));

        Info<< "Selecting laminar thermophysical transport model "
            << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown laminar thermophysical transport model "
                << modelType << nl << nl
                << "Available models:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<laminarThermophysicalTransportModel>
        (
            cstrIter()(momentumTransport, thermo)
        );
    }
    else
    {
        Info<< "Selecting default laminar thermophysical transport model "
            << laminarThermophysicalTransportModels::unityLewisFourier<
               BasicThermophysicalTransportModel>::typeName << endl;

        return autoPtr<laminarThermophysicalTransportModel>
        (
            new laminarThermophysicalTransportModels::unityLewisFourier
            <
                BasicThermophysicalTransportModel
            >(momentumTransport, thermo)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicThermophysicalTransportModel>
bool Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::read()
{
    if (BasicThermophysicalTransportModel::read())
    {
        laminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= laminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicThermophysicalTransportModel>
void Foam::laminarThermophysicalTransportModel
<
    BasicThermophysicalTransportModel
>::correct()
{
    BasicThermophysicalTransportModel::correct();
}


// ************************************************************************* //
