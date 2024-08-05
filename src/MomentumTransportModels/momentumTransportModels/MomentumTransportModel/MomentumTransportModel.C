#include "MomentumTransportModelSB.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    class Alpha,
    class Rho,
    class BasicMomentumTransportModel,
    class TransportModel
>
Foam::MomentumTransportModel
<
    Alpha,
    Rho,
    BasicMomentumTransportModel,
    TransportModel
>::MomentumTransportModel
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    BasicMomentumTransportModel
    (
        rho,
        U,
        alphaRhoPhi,
        phi
    ),
    alpha_(alpha),
    transport_(transport)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template
<
    class Alpha,
    class Rho,
    class BasicMomentumTransportModel,
    class TransportModel
>
Foam::autoPtr
<
    Foam::MomentumTransportModel
    <
        Alpha,
        Rho,
        BasicMomentumTransportModel,
        TransportModel
    >
>
Foam::MomentumTransportModel
<
    Alpha,
    Rho,
    BasicMomentumTransportModel,
    TransportModel
>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
{
    const word modelType
    (
        IOdictionary
        (
            momentumTransportModel::readModelDict
            (
                U.db(),
                alphaRhoPhi.group()
            )
        ).lookup("simulationType")
    );

    Info<< "Selecting turbulence model type " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown MomentumTransportModel type "
            << modelType << nl << nl
            << "Valid MomentumTransportModel types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<MomentumTransportModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport)
    );
}
