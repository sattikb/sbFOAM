// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Thermo>
inline Foam::viscSB<Thermo>::viscSB
(
    const Thermo& t,
    const scalar muInf,
    const scalar mu0,
    const scalar kSB,
    const scalar nSB,
    const scalar kappa
)
:
    Thermo(t),
    muInf_(muInf),
    mu0_(mu0),
    kSB_(kSB),
    nSB_(nSB),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::viscSB<Thermo>::viscSB
(
    const word& name,
    const viscSB& vsb
)
:
    Thermo(name, vsb),
    muInf_(vsb.muInf_),
    mu0_(vsb.mu0_),
    kSB_(vsb.kSB_),
    nSB_(vsb.nSB_),
    kappa_(vsb.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::viscSB<Thermo>>
Foam::viscSB<Thermo>::clone() const
{
    return autoPtr<viscSB<Thermo>>
    (
        new viscSB<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::viscSB<Thermo>>
Foam::viscSB<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<viscSB<Thermo>>
    (
        new viscSB<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::viscSB<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    Info<<"DO NOT USE THIS. ERR ERR ERR ERR ERR SATTIK READ THIS"<<endl;
    return muInf_;
}

template<class Thermo>
inline Foam::scalar Foam::viscSB<Thermo>::mu3ArgSB
(
    const scalar p,
    const scalar T,
    const scalar sr
) const
{
    
    return muInf_ + (mu0_ - muInf_) / (1.0 + pow(kSB_*sr, nSB_) ) ;
}

template<class Thermo>
inline Foam::scalar Foam::viscSB<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::viscSB<Thermo>::operator+=
(
    const viscSB<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        muInf_ = Y1*muInf_ + Y2*st.muInf_;
        mu0_ = Y1*mu0_ + Y2*st.mu0_;
        kSB_ = Y1*kSB_ + Y2*st.kSB_;
        nSB_ = Y1*nSB_ + Y2*st.nSB_;
        kappa_ = Y1*kappa_ + Y2*st.kappa_;
    }
}


template<class Thermo>
inline void Foam::viscSB<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::viscSB<Thermo> Foam::operator+
(
    const viscSB<Thermo>& vsb1,
    const viscSB<Thermo>& vsb2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(vsb1) + static_cast<const Thermo&>(vsb2)
    );

    if (mag(t.Y()) < small)
    {
        return viscSB<Thermo>
        (
            t,
            0,
            vsb1.muInf_,
	    vsb1.mu0_,
	    vsb1.kSB_,
	    vsb1.nSB_,
	    vsb1.kappa_ 
        );
    }
    else
    {
        scalar Y1 = vsb1.Y()/t.Y();
        scalar Y2 = vsb2.Y()/t.Y();

        return viscSB<Thermo>
        (
            t,
            Y1*vsb1.muInf_ + Y2*vsb2.muInf_,
            Y1*vsb1.mu0_ + Y2*vsb2.mu0_,
            Y1*vsb1.kSB_ + Y2*vsb2.kSB_,
            Y1*vsb1.nSB_ + Y2*vsb2.nSB_,
            Y1*vsb1.kappa_ + Y2*vsb2.kappa_
        );
    }
}


template<class Thermo>
inline Foam::viscSB<Thermo> Foam::operator*
(
    const scalar s,
    const viscSB<Thermo>& vsb
)
{
    return viscSB<Thermo>
    (
        s*static_cast<const Thermo&>(vsb),
        vsb.muInf_,
	vsb.mu0_,
	vsb.kSB_,
	vsb.nSB_,
        vsb.kappa_
    );
}
