#ifndef sbRhoEOS_H
#define sbRhoEOS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class sbRhoEOS;

template<class Specie>
inline sbRhoEOS<Specie> operator+
(
    const sbRhoEOS<Specie>&,
    const sbRhoEOS<Specie>&
);

template<class Specie>
inline sbRhoEOS<Specie> operator*
(
    const scalar,
    const sbRhoEOS<Specie>&
);

template<class Specie>
inline sbRhoEOS<Specie> operator==
(
    const sbRhoEOS<Specie>&,
    const sbRhoEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const sbRhoEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                          Class sbRhoEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class sbRhoEOS
:
    public Specie
{
    // Private Data

        scalar rhoRef_;
        scalar rhoMax_;
        scalar pRef_;
        scalar xi_;


public:

    // Constructors

        //- Construct from components
        inline sbRhoEOS
	(
	 	const Specie& sp, 
		const scalar rhoRef,
		const scalar rhoMax,
		const scalar pRef,
		const scalar xi
	);

        //- Construct from dictionary
        sbRhoEOS(const dictionary& dict);

        //- Construct as named copy
        inline sbRhoEOS(const word& name, const sbRhoEOS&);

        //- Construct and return a clone
        inline autoPtr<sbRhoEOS> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sbRhoEOS<" + word(Specie::typeName_()) + '>';
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

        inline void operator+=(const sbRhoEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend sbRhoEOS operator+ <Specie>
        (
            const sbRhoEOS&,
            const sbRhoEOS&
        );

        friend sbRhoEOS operator* <Specie>
        (
            const scalar s,
            const sbRhoEOS&
        );

        friend sbRhoEOS operator== <Specie>
        (
            const sbRhoEOS&,
            const sbRhoEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const sbRhoEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sbRhoEOSI.H"

#ifdef NoRepository
    #include "sbRhoEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
