// CREATES THE TYPEDEFS FOR TURBULENCE, RAS & LES MODELS 
// FOR COMPRESSIBLE FLOWS USING STANDARD LAMINAR FLOW MODELS

#ifndef dynamicMomentumTransportModel_H
#define dynamicMomentumTransportModel_H

#include "CompressibleMomentumTransportModelSB.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "dynamicTransportModelSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressible
    {
        typedef CompressibleMomentumTransportModel<dynamicTransportModel>
            momentumTransportModel;
        typedef laminarModel<momentumTransportModel> laminarModel;
        typedef RASModel<momentumTransportModel> RASModel;
        typedef LESModel<momentumTransportModel> LESModel;

        template<class BasicCompressibleMomentumTransportModel>
        autoPtr<BasicCompressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleMomentumTransportModel::
                transportModel& transport
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicMomentumTransportModel.C"
#endif

#endif
