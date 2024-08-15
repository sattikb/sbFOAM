//RETURNS THE cP/cV (DEPENDS ON THE thermo IN thermophysicalProperties),
// SENSIBLE ENTHALPY, AND THE TEMPERATURE FROM ENTHALPY

#ifndef sensibleEnthalpySB_H
#define sensibleEnthalpySB_H

namespace Foam
{

template<class Thermo>
class sensibleEnthalpy
{

public:

    // Constructors

        //- Construct
        sensibleEnthalpy()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleEnthalpy";
        }

        // Fundamental properties

            static bool enthalpy()
            {
                return true;
            }

            static word energyName()
            {
                return "h";
            }

            // Heat capacity at constant pressure [J/kg/K]
            scalar Cpv
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Cp(p, T);
            }

            // Sensible enthalpy [J/kg]
            scalar HE
            (
                const Thermo& thermo,
                const scalar p,
                const scalar T
            ) const
            {
                return thermo.Hs(p, T);
            }

            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar THE
            (
                const Thermo& thermo,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
		//THs IS DEFINED IN thermoI.H
                return thermo.THs(h, p, T0);
            }
};


} // End namespace Foam

#endif
