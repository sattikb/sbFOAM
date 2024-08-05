// THIS CLASS REPRESENT THE FIELD OF 1, TO AVOID UNNEXESSARY MANIPULATION 
// FOR OBJECTS WHICH ARE 0.

#ifndef oneFieldSB_H
#define oneFieldSB_H

#include "oneSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class oneField Declaration
\*---------------------------------------------------------------------------*/

class oneField
:
    public one
{

public:

    // Public Typedefs

        typedef oneField FieldType;


    // Constructors

        //- Construct null
        oneField()
        {}


    // Member Operators

        inline one operator[](const label) const;

        inline oneField field() const;
};


inline const oneField& operator*(const oneField&, const oneField&);

inline const oneField& operator/(const oneField&, const oneField&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "oneFieldI.H"

#endif
