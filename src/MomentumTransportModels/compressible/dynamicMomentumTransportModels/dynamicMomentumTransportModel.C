#include "dynamicMomentumTransportModelSB.H"

namespace Foam
{
    namespace compressible
    {
        template<class BasicCompressibleMomentumTransportModel>
        autoPtr<BasicCompressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleMomentumTransportModel::
                transportModel& transport
        )
        {
            return BasicCompressibleMomentumTransportModel::New
            (
                geometricOneField(),
                rho,
                U,
                phi,
                phi,
                transport
            );
        }
    }
}
