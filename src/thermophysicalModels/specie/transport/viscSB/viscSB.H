// TRANSPORT OF MU AND KAPPA FOR IFG
//    transport
//    {
//        mu		1.82e-05;
//        kappa		0.71;
//    }

#ifndef viscSB_H
#define viscSB_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class viscSB;

template<class Thermo>
inline viscSB<Thermo> operator+
(
    const viscSB<Thermo>&,
    const viscSB<Thermo>&
);

template<class Thermo>
inline viscSB<Thermo> operator*
(
    const scalar,
    const viscSB<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const viscSB<Thermo>&
);


/*---------------------------------------------------------------------------*\
                       Class viscSB Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class viscSB
:
    public Thermo
{
    // Private Data

        //- Constant dynamic viscosity at Infinite Shear Stress [Pa.s]
        scalar muInf_;
        //- Constant dynamic viscosity at Zero Shear Stress [Pa.s]
        scalar mu0_;
        //- Constant k
        scalar kSB_;
        //- Constant n
        scalar nSB_;
        //- Constant thermal conductivity [W/mK]
        scalar kappa_;


    // Private Member Functions

        //- Construct from components
        inline viscSB
        (
            const Thermo& t,
            const scalar muInf,
            const scalar mu0,
            const scalar kSB,
            const scalar nSB,
            const scalar kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline viscSB(const word&, const viscSB&);

        //- Construct from dictionary
        viscSB(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<viscSB> clone() const;

        // Selector from dictionary
        inline static autoPtr<viscSB> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "viscSB<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;
	inline scalar mu3ArgSB(const scalar p, const scalar T, const scalar sr) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const viscSB&);
        inline void operator*=(const scalar);


    // Friend operators

        friend viscSB operator+ <Thermo>
        (
            const viscSB&,
            const viscSB&
        );

        friend viscSB operator* <Thermo>
        (
            const scalar,
            const viscSB&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const viscSB&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "viscSBI.H"

#ifdef NoRepository
    #include "viscSB.C"
#endif

#endif
