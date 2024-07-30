// TEMPLATE OF IDENTITY AND DUAL SPACE IDENTITY TENSORS DERIVED FROM 
// SPHERICALTENSOR

#ifndef IdentitySB_H
#define IdentitySB_H

#include "SphericalTensorSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Identity Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Identity
:
    public SphericalTensor<Cmpt>
{

public:

    //- Construct initialising the SphericalTensor to 1
    Identity()
    :
        SphericalTensor<Cmpt>(1)
    {}

    //- The identity type in the dual space
    class dual
    :
        public SphericalTensor<Cmpt>
    {
        //- Construct initialising the SphericalTensor to 1
        dual()
        :
            SphericalTensor<Cmpt>(1)
        {}
    };

    //- Return the identity in the dual space
    inline dual operator*()
    {
        return dual();
    }

    //- Return 1 for label
    inline explicit operator label() const
    {
        return 1;
    }

    //- Return 1 for scalar
    inline explicit operator scalar() const
    {
        return 1;
    }
};


// Global Identity tensor
static const Identity<scalar> I;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
