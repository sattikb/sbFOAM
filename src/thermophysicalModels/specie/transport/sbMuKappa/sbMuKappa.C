#include "sbMuKappa.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Thermo>
Foam::sbMuKappa<Thermo>::sbMuKappa(const dictionary& dict)
:
    Thermo(dict),
    mu_(dict.subDict("transport").lookup<scalar>("mu")),
    kappa_(dict.subDict("transport").lookup<scalar>("kappa"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Thermo>
void Foam::sbMuKappa<Thermo>::sbMuKappa::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("mu", mu_);
    dict.add("kappa", kappa_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const sbMuKappa<Thermo>& sbmk)
{
    sbmk.write(os);
    return os;
}

