#include "ifgRhoEOS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::ifgRhoEOS<Specie>::ifgRhoEOS
(
    const Specie& sp,
    const scalar MC
)
:
    Specie(sp),
    MC_(MC)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::ifgRhoEOS<Specie>::ifgRhoEOS
(
    const word& name,
    const ifgRhoEOS<Specie>& ifgrt
)
:
    Specie(name, ifgrt),
    MC_(ifgrt.MC_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::ifgRhoEOS<Specie>>
Foam::ifgRhoEOS<Specie>::clone() const
{
    return autoPtr<ifgRhoEOS<Specie>>(new ifgRhoEOS<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::gammaSB
(
    const scalar x,
    const scalar y,
    const scalar kk
) const
{
    scalar g = 0;
    scalar d = x;
    for (scalar k = 0; k <= kk; ++k)
    {
        g = g + pow(y, x) * exp(-y) * pow(y, k) / d;
        d = d * (x + k + scalar(1.));
    }
    return g;
}

template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::rho(scalar p, scalar T) const
{
    const scalar a = -9.430631e-14;
    const scalar b =  2.804769e-08;
    const scalar c = -2.522709e-03;
    const scalar d =  258.879724;

    const scalar e = 7.4e-3;
    const scalar f = -0.826968;

    const scalar ga = 0.3796;
    const scalar gb = 9.622e-8;

    const scalar mcFact = scalar(1.0) - MC_;

    if (p<scalar(-32375.43))			//68949.57
    {
	const scalar pp = -32375.43;
	const scalar rhoSum = a * pow(pp,3) + b * pow(pp,2) + c * pp + d;
    	return rhoSum / mcFact;
    }
    else if (p>scalar(-32375.43) && p<scalar(25.0))	//68949.57, 101350
    {
	const scalar rhoSum = a * pow(p,3) + b * pow(p,2) + c * p + d;
    	return rhoSum / mcFact;
    }
    else
    {
	const scalar rhoSum = scalar(1.0) / (e * exp( f * gammaSB(ga, gb*p, 10) ) );
	return rhoSum / mcFact;
    }


}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::Cp(scalar p, scalar T) const
{
    return 0;	//WILL BE TAKEN FROM ~hConst
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::E(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::Cv(scalar p, scalar T) const
{
    return 0;	//WILL BE TAKEN FROM ~eConst
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::Sp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::Sv(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::psi(scalar p, scalar T) const
{
    //This is drho/dp
    const scalar a = -9.430631e-14;
    const scalar b =  2.804769e-08;
    const scalar c = -2.522709e-03;

    const scalar e = 7.4e-3;
    const scalar f = -0.826968;

    const scalar ga = 0.3796;
    const scalar gb = 9.622e-8;

    const scalar mcFact = scalar(1.0) - MC_;

    if (p<scalar(-32375.43))			//68949.57
    {
	const scalar pp = -32375.43;		//68949.57
	const scalar drdpSum = 3*a*pow(pp,2) + 2*b*pp + c;
    	return drdpSum / mcFact;
    }
    else if (p>scalar(-32375.43) && p<scalar(25.0))	//68949.57, 101350
    {
	const scalar drdpSum = 3*a*pow(p,2) + 2*b*p + c;
    	return drdpSum / mcFact;
    }
    else
    {
	const scalar rhoSum = scalar(1.0) / (e * exp( f * gammaSB(ga, gb*p, 10) ) );
	return - rhoSum / mcFact * f * gb * pow (gb*p, e-scalar(1.0)) * exp(-gb*p);
    }

}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::Z(scalar p, scalar T) const
{
    const scalar a = -9.430631e-14;
    const scalar b =  2.804769e-08;
    const scalar c = -2.522709e-03;
    const scalar d =  258.879724;

    const scalar e = 7.4e-3;
    const scalar f = -0.826968;

    const scalar ga = 0.3796;
    const scalar gb = 9.622e-8;

    const scalar mcFact = scalar(1.0) - MC_;

    if (p<scalar(-32375.43))			//68949.57
    {
	const scalar pp = -32375.43;		//68949.57
	const scalar rhoSum = a * pow(pp,3) + b * pow(pp,2) + c * pp + d;
    	const scalar thisrho = rhoSum / mcFact;
	return p/(thisrho*this->R()*T);
    }
    else if (p>scalar(-32375.43) && p<scalar(25.0))	//68949.57, 101350
    {
	const scalar rhoSum = a * pow(p,3) + b * pow(p,2) + c * p + d;
    	const scalar thisrho = rhoSum / mcFact;
	return p/(thisrho*this->R()*T);
    }
    else
    {
	const scalar rhoSum = scalar(1.0) / (e * exp( f * gammaSB(ga, gb*p, 10) ) );
    	const scalar thisrho = rhoSum / mcFact;
	return p/(thisrho*this->R()*T);
    }
}


template<class Specie>
inline Foam::scalar Foam::ifgRhoEOS<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::ifgRhoEOS<Specie>::operator+=(const ifgRhoEOS<Specie>& ifgrt)
{
    const scalar Y1 = this->Y();
    Specie::operator+=(ifgrt);

    if (mag(this->Y()) > small)
    {
        MC_ = this->Y()/(Y1/MC_ + ifgrt.Y()/ifgrt.MC_);
    }
}


template<class Specie>
inline void Foam::ifgRhoEOS<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::ifgRhoEOS<Specie> Foam::operator+
(
    const ifgRhoEOS<Specie>& ifgrt1,
    const ifgRhoEOS<Specie>& ifgrt2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ifgrt1)
      + static_cast<const Specie&>(ifgrt2)
    );

    if (mag(sp.Y()) < small)
    {
        return ifgRhoEOS<Specie>
        (
            sp,
            ifgrt1.MC_
        );
    }
    else
    {
        return ifgRhoEOS<Specie>
        (
            sp,
            sp.Y()/(ifgrt1.Y()/ifgrt1.MC_ + ifgrt2.Y()/ifgrt2.MC_)
        );
    }
}


template<class Specie>
inline Foam::ifgRhoEOS<Specie> Foam::operator*
(
    const scalar s,
    const ifgRhoEOS<Specie>& ifgrt
)
{
    return ifgRhoEOS<Specie>
    (
	s*static_cast<const Specie&>(ifgrt), 
	ifgrt.MC_
    );
}


template<class Specie>
inline Foam::ifgRhoEOS<Specie> Foam::operator==
(
    const ifgRhoEOS<Specie>& ifgrt1,
    const ifgRhoEOS<Specie>& ifgrt2
)
{
    return ifgRhoEOS<Specie>
    (
        //Might need some modifications
        static_cast<const Specie&>(ifgrt1) == static_cast<const Specie&>(ifgrt2),
        NaN
    );
}
