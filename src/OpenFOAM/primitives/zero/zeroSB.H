// CLASS REPRESENTING THE CONCEPT OF 0. THIS AVOIDS MANIPULATIONS ON FIELDS
// KNOWN TO BE 0 AT COMPILE-TIME.

#ifndef zeroSB_H
#define zeroSB_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class zero Declaration
\*---------------------------------------------------------------------------*/

class zero
{
public:

    typedef zero value_type;

    // Constructors

        //- Construct null
        zero()
        {}


    // Member Operators

        //- Return 0 for bool
        inline operator bool() const
        {
            return 0;
        }

        //- Return 0 for label
        inline operator label() const
        {
            return 0;
        }

        //- Return 0 for float
        inline operator float() const
        {
            return 0;
        }

        //- Return 0 for double
        inline operator double() const
        {
            return 0;
        }

        //- Return 0 for double
        inline operator long double() const
        {
            return 0;
        }
};


// Global zero
static const zero Zero;


} // End namespace Foam

#include "zeroI.H"

#endif
