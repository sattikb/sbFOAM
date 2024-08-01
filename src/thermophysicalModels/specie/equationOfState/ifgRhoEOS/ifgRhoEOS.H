#ifndef ifgRhoEOS_H
#define ifgRhoEOS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class ifgRhoEOS;

template<class Specie>
inline ifgRhoEOS<Specie> operator+
(
    const ifgRhoEOS<Specie>&,
    const ifgRhoEOS<Specie>&
);

template<class Specie>
inline ifgRhoEOS<Specie> operator*
(
    const scalar,
    const ifgRhoEOS<Specie>&
);

template<class Specie>
inline ifgRhoEOS<Specie> operator==
(
    const ifgRhoEOS<Specie>&,
    const ifgRhoEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const ifgRhoEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                          Class ifgRhoEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class ifgRhoEOS
:
    public Specie
{
    // Private Data

        scalar MC_;

	scalar gammaSB(scalar x, scalar y, scalar kk) const;


public:

    // Constructors

        //- Construct from components
        inline ifgRhoEOS
	(
	 	const Specie& sp, 
		const scalar MC
	);

        //- Construct from dictionary
        ifgRhoEOS(const dictionary& dict);

        //- Construct as named copy
        inline ifgRhoEOS(const word& name, const ifgRhoEOS&);

        //- Construct and return a clone
        inline autoPtr<ifgRhoEOS> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "ifgRhoEOS<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            inline scalar Sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            inline scalar Sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const ifgRhoEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend ifgRhoEOS operator+ <Specie>
        (
            const ifgRhoEOS&,
            const ifgRhoEOS&
        );

        friend ifgRhoEOS operator* <Specie>
        (
            const scalar s,
            const ifgRhoEOS&
        );

        friend ifgRhoEOS operator== <Specie>
        (
            const ifgRhoEOS&,
            const ifgRhoEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const ifgRhoEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ifgRhoEOSI.H"

#ifdef NoRepository
    #include "ifgRhoEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
