// FROM Vector, A SCALAR TYPEDEF IS CREATED CALLED vector.

#ifndef vector_H
#define vector_H

#include "VectorSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Vector<scalar> vector;


//- Data associated with vector type are contiguous
template<>
inline bool contiguous<vector>() {return true;}


template<class Type>
class flux
:
    public innerProduct<vector, Type>
{};


template<>
class flux<scalar>
{
public:

    typedef scalar type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
