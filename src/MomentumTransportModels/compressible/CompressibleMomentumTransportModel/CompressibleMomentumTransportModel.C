#include "CompressibleMomentumTransportModelSB.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::CompressibleMomentumTransportModel<TransportModel>::
CompressibleMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport
)
:
    MomentumTransportModel
    <
        geometricOneField,
        volScalarField,
        compressibleMomentumTransportModel,
        transportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::CompressibleMomentumTransportModel<TransportModel>>
Foam::CompressibleMomentumTransportModel<TransportModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport
)
{
    return autoPtr<CompressibleMomentumTransportModel>
    (
        static_cast<CompressibleMomentumTransportModel*>(
        MomentumTransportModel
        <
            geometricOneField,
            volScalarField,
            compressibleMomentumTransportModel,
            transportModel
        >::New
        (
            geometricOneField(),
            rho,
            U,
            phi,
            phi,
            transport
        ).ptr())
    );
}
