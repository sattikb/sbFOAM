#include "viscSB.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Thermo>
Foam::viscSB<Thermo>::viscSB(const dictionary& dict)
:
    Thermo(dict),
    muInf_(dict.subDict("transport").lookup<scalar>("muInf")),
    mu0_(dict.subDict("transport").lookup<scalar>("mu0")),
    kSB_(dict.subDict("transport").lookup<scalar>("k")),
    nSB_(dict.subDict("transport").lookup<scalar>("n")),
    kappa_(dict.subDict("transport").lookup<scalar>("kappa"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::viscSB<Thermo>::viscSB::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("muInf", muInf_);
    dict.add("mu0", mu0_);
    dict.add("kSB", kSB_);
    dict.add("nSB", nSB_);
    dict.add("kappa", kappa_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const viscSB<Thermo>& vsb)
{
    vsb.write(os);
    return os;
}

