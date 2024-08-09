// TYPEDEFS FOR LAMINAR AND TURBULENT THERMOPHYSICAL TRANSPORT MODELS FOR 
// FLUIDTHERMO FLOW.
//  Foam::compressible::thermophysicalTransportModel
//  Foam::compressible::laminarThermophysicalTransportModel
//  Foam::compressible::RASThermophysicalTransportModel
//  Foam::compressible::LESThermophysicalTransportModel

#ifndef fluidThermophysicalTransportModelSB_H
#define fluidThermophysicalTransportModelSB_H

//#include "fluidThermo.H"
#include "laminarThermophysicalTransportModelSB.H"
#include "RASThermophysicalTransportModel.H"
#include "LESThermophysicalTransportModel.H"
#include "dynamicMomentumTransportModelSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ThermophysicalTransportModel
    <
        compressible::momentumTransportModel,
        fluidThermo
    > fluidThermophysicalTransportModel;

    typedef laminarThermophysicalTransportModel
    <
        fluidThermophysicalTransportModel
    > laminarFluidThermophysicalTransportModel;

    typedef RASThermophysicalTransportModel
    <
        fluidThermophysicalTransportModel
    > RASFluidThermophysicalTransportModel;

    typedef LESThermophysicalTransportModel
    <
        fluidThermophysicalTransportModel
    > LESFluidThermophysicalTransportModel;
}

#endif
