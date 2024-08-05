#include "oneFieldFieldSB.H"

inline Foam::oneField Foam::oneFieldField::operator[](const label) const
{
    return oneField();
}
