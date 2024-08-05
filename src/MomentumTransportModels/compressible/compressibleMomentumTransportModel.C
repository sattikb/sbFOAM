#include "compressibleMomentumTransportModelSB.H"
#include "surfaceInterpolate.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleMomentumTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleMomentumTransportModel::compressibleMomentumTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi
)
:
    momentumTransportModel
    (
        U,
        alphaRhoPhi,
        phi
    ),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::compressibleMomentumTransportModel::phi() const
{
    if (phi_.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        return phi_;	//volumetric flow rate m3/s
    }
    else
    {
        return phi_/fvc::interpolate(rho_);
    }
}
