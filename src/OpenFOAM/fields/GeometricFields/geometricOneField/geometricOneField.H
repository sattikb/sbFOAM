// THIS CLASS ENSURES ANY OBJECT THAT WILL BE 1 (SCALAR, FIELD, FIELDFIELD) 
// ALL REMAIN ONE TO AVOID UNNECESSARY MANIPULATIONS.
//  Used for example as the density argument to a function written for
//  compressible to be used for incompressible flow.

#ifndef geometricOneField_H
#define geometricOneField_H

#include "oneFieldFieldSB.H"
#include "dimensionSetSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class geometricOneField Declaration
\*---------------------------------------------------------------------------*/

class geometricOneField
:
    public one
{

public:

    // Public Typedefs

        typedef oneField Internal;
        typedef oneField Patch;
        typedef oneFieldField Boundary;
        typedef one cmptType;

    // Constructors

        //- Construct null
        geometricOneField()
        {}


    // Member Operators

        inline const dimensionSet& dimensions() const;

        inline one operator[](const label) const;

        inline oneField field() const;

        inline oneField oldTime() const;

        inline Internal operator()() const;

        inline Internal v() const;

        inline typename Internal::FieldType primitiveField() const;

        inline Boundary boundaryField() const;
};


inline const geometricOneField& operator*
(
    const geometricOneField&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator*
(
    const geometricOneField::Internal&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator*
(
    const geometricOneField&,
    const geometricOneField::Internal&
);

inline const geometricOneField& operator/
(
    const geometricOneField&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator/
(
    const geometricOneField::Internal&,
    const geometricOneField&
);

inline const geometricOneField::Internal& operator/
(
    const geometricOneField&,
    const geometricOneField::Internal&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "geometricOneFieldI.H"

#endif
