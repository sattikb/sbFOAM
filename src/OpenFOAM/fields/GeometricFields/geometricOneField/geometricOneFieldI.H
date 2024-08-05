// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::geometricOneField::dimensions() const
{
    return dimless;
}


inline Foam::one Foam::geometricOneField::operator[](const label) const
{
    return one();
}


inline Foam::oneField Foam::geometricOneField::field() const
{
    return oneField();
}


inline Foam::oneField Foam::geometricOneField::oldTime() const
{
    return oneField();
}


inline Foam::geometricOneField::Internal
Foam::geometricOneField::operator()() const
{
    return Internal();
}


inline Foam::geometricOneField::Internal
Foam::geometricOneField::v() const
{
    return Internal();
}


inline typename Foam::geometricOneField::Internal::FieldType
Foam::geometricOneField::primitiveField() const
{
    return typename Internal::FieldType();
}


inline Foam::geometricOneField::Boundary
Foam::geometricOneField::boundaryField() const
{
    return Boundary();
}


inline const Foam::geometricOneField& Foam::operator*
(
    const geometricOneField& gof,
    const geometricOneField&
)
{
    return gof;
}


inline const Foam::geometricOneField::Internal& Foam::operator*
(
    const geometricOneField::Internal& of,
    const geometricOneField&
)
{
    return of;
}


inline const Foam::geometricOneField::Internal& Foam::operator*
(
    const geometricOneField&,
    const geometricOneField::Internal& of
)
{
    return of;
}


inline const Foam::geometricOneField& Foam::operator/
(
    const geometricOneField& gof,
    const geometricOneField&
)
{
    return gof;
}


inline const Foam::geometricOneField::Internal& Foam::operator/
(
    const geometricOneField::Internal& of,
    const geometricOneField&
)
{
    return of;
}


inline const Foam::geometricOneField::Internal& Foam::operator/
(
    const geometricOneField&,
    const geometricOneField::Internal& of
)
{
    return of;
}


// ************************************************************************* //
