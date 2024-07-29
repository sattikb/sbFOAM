#include "sbRhoEOS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::sbRhoEOS<Specie>::sbRhoEOS
(
    const Specie& sp,
    const scalar rhoRef,
    const scalar rhoMax,
    const scalar pRef,
    const scalar xi
)
:
    Specie(sp),
    rhoRef_(rhoRef),
    rhoMax_(rhoMax),
    pRef_(pRef),
    xi_(xi)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::sbRhoEOS<Specie>::sbRhoEOS
(
    const word& name,
    const sbRhoEOS<Specie>& sbrt
)
:
    Specie(name, sbrt),
    rhoRef_(sbrt.rhoRef_),
    rhoMax_(sbrt.rhoMax_),
    pRef_(sbrt.pRef_),
    xi_(sbrt.xi_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::sbRhoEOS<Specie>>
Foam::sbRhoEOS<Specie>::clone() const
{
    return autoPtr<sbRhoEOS<Specie>>(new sbRhoEOS<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::rho(scalar p, scalar T) const
{
    // Duncan et al - Pressure driven flow of Lignocellulosic Biomass
    // A compressible Bingham fluid (2018)
    const scalar powSB = scalar(1.0) - pow( (p/pRef_), (scalar(1.0)-xi_) );
    return 100;//rhoRef_*pow( (rhoMax_/rhoRef_), powSB );
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::H(scalar p, scalar T) const
{
    //I modified this as well
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::Cv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::Sp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::Sv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::psi(scalar p, scalar T) const
{
    // Need to check this
    const scalar powSB  = scalar(1.0) - pow( (p/pRef_), (scalar(1.0)-xi_) );
    const scalar rhoSB  =  rhoRef_*pow( (rhoMax_/rhoRef_), powSB );
    const scalar fact   = (xi_-scalar(1.0))*pow( (p/pRef_), -xi_);
    const scalar compSB = (rhoSB/pRef_)*log(rhoMax_/rhoRef_)*fact;
    return compSB;
//    return 0;
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::Z(scalar p, scalar T) const
{
    // Need to check this
    const scalar powSB = scalar(1.0) - pow( (p/pRef_), (scalar(1.0)-xi_) );
    const scalar rhoSB =  rhoRef_*pow( (rhoMax_/rhoRef_), powSB );
    return p/(rhoSB*this->R()*T);
//    return p/(rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::sbRhoEOS<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::sbRhoEOS<Specie>::operator+=(const sbRhoEOS<Specie>& sbrt)
{
    const scalar Y1 = this->Y();
    Specie::operator+=(sbrt);

    if (mag(this->Y()) > small)
    {
        rhoRef_ = this->Y()/(Y1/rhoRef_ + sbrt.Y()/sbrt.rhoRef_);
        rhoMax_ = this->Y()/(Y1/rhoMax_ + sbrt.Y()/sbrt.rhoMax_);
        pRef_ = this->Y()/(Y1/pRef_ + sbrt.Y()/sbrt.pRef_);
        xi_ = this->Y()/(Y1/xi_ + sbrt.Y()/sbrt.xi_);
    }
}


template<class Specie>
inline void Foam::sbRhoEOS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::sbRhoEOS<Specie> Foam::operator+
(
    const sbRhoEOS<Specie>& sbrt1,
    const sbRhoEOS<Specie>& sbrt2
)
{
    Specie sp
    (
        static_cast<const Specie&>(sbrt1)
      + static_cast<const Specie&>(sbrt2)
    );

    if (mag(sp.Y()) < small)
    {
        return sbRhoEOS<Specie>
        (
            sp,
            sbrt1.rhoRef_,
            sbrt1.rhoMax_,
            sbrt1.pRef_,
            sbrt1.xi_
        );
    }
    else
    {
        return sbRhoEOS<Specie>
        (
            sp,
            sp.Y()/(sbrt1.Y()/sbrt1.rhoRef_ + sbrt2.Y()/sbrt2.rhoRef_),
            sp.Y()/(sbrt1.Y()/sbrt1.rhoMax_ + sbrt2.Y()/sbrt2.rhoMax_),
            sp.Y()/(sbrt1.Y()/sbrt1.pRef_ + sbrt2.Y()/sbrt2.pRef_),
            sp.Y()/(sbrt1.Y()/sbrt1.xi_ + sbrt2.Y()/sbrt2.xi_)
        );
    }
}


template<class Specie>
inline Foam::sbRhoEOS<Specie> Foam::operator*
(
    const scalar s,
    const sbRhoEOS<Specie>& sbrt
)
{
    return sbRhoEOS<Specie>
    (
	s*static_cast<const Specie&>(sbrt), 
	sbrt.rhoRef_,
	sbrt.rhoMax_,
	sbrt.pRef_,
	sbrt.xi_
    );
}


template<class Specie>
inline Foam::sbRhoEOS<Specie> Foam::operator==
(
    const sbRhoEOS<Specie>& sbrt1,
    const sbRhoEOS<Specie>& sbrt2
)
{
    return sbRhoEOS<Specie>
    (
        //Might need some modifications
        static_cast<const Specie&>(sbrt1) == static_cast<const Specie&>(sbrt2),
        NaN,
	NaN,
	NaN,
	NaN
    );
}
