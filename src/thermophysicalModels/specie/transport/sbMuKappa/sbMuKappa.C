#include "sbMuKappa.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Thermo>
Foam::constTransport<Thermo>::constTransport(const dictionary& dict)
:
    Thermo(dict),
    mu_(dict.subDict("transport").lookup<scalar>("mu")),
    rPr_(1.0/dict.subDict("transport").lookup<scalar>("Pr"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::constTransport<Thermo>::constTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("mu", mu_);
    dict.add("Pr", 1.0/rPr_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const constTransport<Thermo>& sbmk)
{
    sbmk.write(os);
    return os;
}

