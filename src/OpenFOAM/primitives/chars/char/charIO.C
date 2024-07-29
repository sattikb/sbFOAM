// DEFINING THE INSTANCES FROM charSB.H WITH FURTHER FWD DECLARED METHODS

#include "charSB.H"
#include "IOstreams.H"

char Foam::readChar(Istream& is)
{
   char c;
   is.read(c);
   return c;
}


Foam::Istream& Foam::operator>>(Istream& is, char& c)
{
    is.read(c);
    is.check("Istream& operator>>(Istream&, char&)");
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const char c)
{
    os.write(c);
    os.check("Ostream& operator<<(Ostream&, const char)");
    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const char* s)
{
    os.write(s);
    os.check("Ostream& operator<<(Ostream&, const char*)");
    return os;
}
