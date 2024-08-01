#include "ifgRhoEOS.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::ifgRhoEOS<Specie>::ifgRhoEOS(const dictionary& dict)
:
    Specie(dict),
    MC_(dict.subDict("equationOfState").lookup<scalar>("MC"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::ifgRhoEOS<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("MC", MC_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const ifgRhoEOS<Specie>& ifgrt)
{
    ifgrt.write(os);
    return os;
}
