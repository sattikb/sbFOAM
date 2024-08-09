// TRANSPORT OF MU AND KAPPA FOR IFG
//    transport
//    {
//        mu		1.82e-05;
//        kappa		0.71;
//    }

#ifndef sbMuKappa_H
#define sbMuKappa_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class sbMuKappa;

template<class Thermo>
inline sbMuKappa<Thermo> operator+
(
    const sbMuKappa<Thermo>&,
    const sbMuKappa<Thermo>&
);

template<class Thermo>
inline sbMuKappa<Thermo> operator*
(
    const scalar,
    const sbMuKappa<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const sbMuKappa<Thermo>&
);


/*---------------------------------------------------------------------------*\
                       Class sbMuKappa Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class sbMuKappa
:
    public Thermo
{
    // Private Data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Constant thermal conductivity [W/mK]
        scalar kappa_;


    // Private Member Functions

        //- Construct from components
        inline sbMuKappa
        (
            const Thermo& t,
            const scalar mu,
            const scalar kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline sbMuKappa(const word&, const sbMuKappa&);

        //- Construct from dictionary
        sbMuKappa(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<sbMuKappa> clone() const;

        // Selector from dictionary
        inline static autoPtr<sbMuKappa> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sbMuKappa<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const sbMuKappa&);
        inline void operator*=(const scalar);


    // Friend operators

        friend sbMuKappa operator+ <Thermo>
        (
            const sbMuKappa&,
            const sbMuKappa&
        );

        friend sbMuKappa operator* <Thermo>
        (
            const scalar,
            const sbMuKappa&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const sbMuKappa&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "sbMuKappaI.H"

#ifdef NoRepository
    #include "sbMuKappa.C"
#endif

#endif
