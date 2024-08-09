#include "ThermophysicalTransportModelSB.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class ThermoModel>
Foam::ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>::
ThermophysicalTransportModel
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    thermophysicalTransportModel(momentumTransport),
    momentumTransport_(momentumTransport),
    thermo_(thermo)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class MomentumTransportModel, class ThermoModel>
Foam::autoPtr
<
    Foam::ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>
>
Foam::ThermophysicalTransportModel<MomentumTransportModel, ThermoModel>::New
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
{
    const word modelType
    (
        momentumTransport.lookup("simulationType")
    );

    Info<< "Selecting thermophysical transport type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermophysical transport type "
            << modelType << nl << nl
            << "Available types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ThermophysicalTransportModel>
    (
        cstrIter()(momentumTransport, thermo)
    );
}


// ************************************************************************* //
