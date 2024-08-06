//CONTAINS THE GENERIC CLASSES OF THE MOMENTUM MODELS, BOTH LAMINAR AND 
//TURBULENT

#include "CompressibleMomentumTransportModelSB.H"
#include "dynamicTransportModelSB.H"
#include "makeMomentumTransportModel.H"

#include "laminarModelSB.H"
#include "RASModel.H"
#include "LESModelSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeMomentumTransportModelTypes
(
    geometricOneField,			// ALPHA (VOL FRACTION)
    volScalarField,			// RHO
    compressibleMomentumTransportModel,
    CompressibleMomentumTransportModel,
    dynamicTransportModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (dynamicTransportModelCompressibleMomentumTransportModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (dynamicTransportModelCompressibleMomentumTransportModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (dynamicTransportModelCompressibleMomentumTransportModel, LES, Type)
