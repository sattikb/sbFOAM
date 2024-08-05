// THIS CLASS REPRESENTS A FIELD OF ONEFIELDS THAT ENSURE THE OBJECTS SUPPOSED 
// TO BE 1 STAY AS 1.

#ifndef oneFieldFieldSB_H
#define oneFieldFieldSB_H

#include "oneFieldSB.H"

namespace Foam
{

class oneFieldField
:
    public one
{

public:

    // Constructors

        //- Construct null
        oneFieldField()
        {}


    // Member Operators

        inline oneField operator[](const label) const;
};


} // End namespace Foam

    #include "oneFieldFieldI.H"

#endif
