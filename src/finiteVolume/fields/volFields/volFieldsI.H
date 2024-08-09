template<class Type, template<class> class PatchField, class GeoMesh>
inline
const typename
Foam::GeometricField<Type, PatchField, GeoMesh>::Internal&
Foam::GeometricField<Type, PatchField, GeoMesh>::
v() const
{
    static_assert(isVolMesh<GeoMesh>::value, "Only valid for volFields");
    return *this;
}
