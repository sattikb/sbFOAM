#include "directionSB.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::direction Foam::readDirection(Istream& is)
{
    direction val;
    is >> val;

    return val;
}


Foam::Istream& Foam::operator>>(Istream& is, direction& d)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        d = direction(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorInFunction(is)
            << "wrong token type - expected direction, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, direction&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const direction d)
{
    os.write(label(d));
    os.check("Ostream& operator<<(Ostream&, const direction)");
    return os;
}


std::ostream& Foam::operator<<(std::ostream& os, const direction d)
{
    os << int(d);
    return os;
}
