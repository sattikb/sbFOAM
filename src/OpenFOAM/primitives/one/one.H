// THIS CLASS REPRESENTS A SCALAR VALUE OF 1, WHICH CAN BE USED TO AVOID
// UNNECESSARY ALGEBRAIC MANIPULATIONS FOR OBJECTS KNOWN TO BE EQUAL TO 1 
// DURING COMPILATION.

#ifndef one_H
#define one_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class one Declaration
\*---------------------------------------------------------------------------*/

class one
{
public:

    typedef one value_type;

    // Constructors

        //- Construct null
        one()
        {}


    // Member Operators

        //- Return 1 for label
        inline operator label() const
        {
            return 1;
        }

        //- Return 1 for float
        inline operator float() const
        {
            return 1;
        }

        //- Return 1 for double
        inline operator double() const
        {
            return 1;
        }

        //- Return 1 for double
        inline operator long double() const
        {
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "oneI.H"

#endif
