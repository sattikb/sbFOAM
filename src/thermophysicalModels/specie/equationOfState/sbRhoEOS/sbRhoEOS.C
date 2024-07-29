#include "sbRhoEOS.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::sbRhoEOS<Specie>::sbRhoEOS(const dictionary& dict)
:
    Specie(dict),
    rhoRef_(dict.subDict("equationOfState").lookup<scalar>("rhoRef")),
    rhoMax_(dict.subDict("equationOfState").lookup<scalar>("rhoMax")),
    pRef_(dict.subDict("equationOfState").lookup<scalar>("pRef")),
    xi_(dict.subDict("equationOfState").lookup<scalar>("xi"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::sbRhoEOS<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("rhoRef", rhoRef_);
    dict.add("rhoMax", rhoMax_);
    dict.add("pRef", pRef_);
    dict.add("xi", xi_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const sbRhoEOS<Specie>& sbrt)
{
    sbrt.write(os);
    return os;
}
