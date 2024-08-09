// TEMPLATED ABC FOR THERMOPHYSICAL TRANSPORT MODELS.

#ifndef ThermophysicalTransportModelSB_H
#define ThermophysicalTransportModelSB_H

#include "thermophysicalTransportModelSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

template<class MomentumTransportModel, class ThermoModel>
class ThermophysicalTransportModel
:
    public thermophysicalTransportModel
{

public:

    typedef MomentumTransportModel momentumTransportModel;
    typedef ThermoModel thermoModel;
    typedef typename momentumTransportModel::alphaField alphaField;


protected:

    // Protected data

        const momentumTransportModel& momentumTransport_;
        const thermoModel& thermo_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            ThermophysicalTransportModel,
            dictionary,
            (
                const momentumTransportModel& momentumTransport,
                const thermoModel& thermo
            ),
            (momentumTransport, thermo)
        );


    // Constructors

        //- Construct
        ThermophysicalTransportModel
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Disallow default bitwise copy construction
        ThermophysicalTransportModel
        (
            const ThermophysicalTransportModel&
        ) = delete;


    // Selectors

        //- Return a reference to the selected thermophysical transport model
        static autoPtr<ThermophysicalTransportModel> New
        (
            const momentumTransportModel& transport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~ThermophysicalTransportModel()
    {}


    // Member Functions

        //- Access function to incompressible transport model
        const momentumTransportModel& momentumTransport() const
        {
            return momentumTransport_;
        }

        //- Access function to incompressible transport model
        virtual const thermoModel& thermo() const
        {
            return thermo_;
        }

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return momentumTransport_.alpha();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ThermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "ThermophysicalTransportModel.C"
#endif

#endif
