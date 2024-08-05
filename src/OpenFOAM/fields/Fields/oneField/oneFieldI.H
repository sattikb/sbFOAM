#include "oneFieldSB.H"

inline Foam::one Foam::oneField::operator[](const label) const
{
    return one();
}


inline Foam::oneField Foam::oneField::field() const
{
    return oneField();
}


inline const Foam::oneField& Foam::operator*
(
    const oneField& of,
    const oneField&
)
{
    return of;
}


inline const Foam::oneField& Foam::operator/
(
    const oneField& of,
    const oneField&
)
{
    return of;
}
