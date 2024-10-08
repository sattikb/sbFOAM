// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sbMuKappa<Thermo>::sbMuKappa
(
    const Thermo& t,
    const scalar mu,
    const scalar kappa
)
:
    Thermo(t),
    mu_(mu),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::sbMuKappa<Thermo>::sbMuKappa
(
    const word& name,
    const sbMuKappa& sbmk
)
:
    Thermo(name, sbmk),
    mu_(sbmk.mu_),
    kappa_(sbmk.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::sbMuKappa<Thermo>>
Foam::sbMuKappa<Thermo>::clone() const
{
    return autoPtr<sbMuKappa<Thermo>>
    (
        new sbMuKappa<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::sbMuKappa<Thermo>>
Foam::sbMuKappa<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sbMuKappa<Thermo>>
    (
        new sbMuKappa<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::sbMuKappa<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::sbMuKappa<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::sbMuKappa<Thermo>::operator+=
(
    const sbMuKappa<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        mu_ = Y1*mu_ + Y2*st.mu_;
        kappa_ = Y1*kappa_ + Y2*st.kappa_;
    }
}


template<class Thermo>
inline void Foam::sbMuKappa<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::sbMuKappa<Thermo> Foam::operator+
(
    const sbMuKappa<Thermo>& sbmk1,
    const sbMuKappa<Thermo>& sbmk2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(sbmk1) + static_cast<const Thermo&>(sbmk2)
    );

    if (mag(t.Y()) < small)
    {
        return sbMuKappa<Thermo>
        (
            t,
            0,
            0
        );
    }
    else
    {
        scalar Y1 = sbmk1.Y()/t.Y();
        scalar Y2 = sbmk2.Y()/t.Y();

        return sbMuKappa<Thermo>
        (
            t,
            Y1*sbmk1.mu_ + Y2*sbmk2.mu_,
            Y1*sbmk1.kappa_ + Y2*sbmk2.kappa_
        );
    }
}


template<class Thermo>
inline Foam::sbMuKappa<Thermo> Foam::operator*
(
    const scalar s,
    const sbMuKappa<Thermo>& sbmk
)
{
    return sbMuKappa<Thermo>
    (
        s*static_cast<const Thermo&>(sbmk),
        sbmk.mu_,
        sbmk.kappa_
    );
}
