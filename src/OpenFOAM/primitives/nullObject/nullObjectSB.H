// CLASS NullObject AND CREATES INSTANCES OF IT

#ifndef nullObjectSB_H
#define nullObjectSB_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class nullObject Declaration
\*---------------------------------------------------------------------------*/

class NullObject
{
    //- Private constructor
    NullObject()
    {}

public:

    //- The unique null object
    static const NullObject nullObject;

    //- Disallow default bitwise copy construction
    NullObject(const NullObject&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const NullObject&) = delete;
};


//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


//- Return reference to the nullObject of type T
template<class T>
inline const T& NullObjectRef();

//- Return pointer to the nullObject of type T
template<class T>
inline const T* NullObjectPtr();

//- Return reference to the nullObject of type T
template<class T>
inline T&& NullObjectMove();


//- Return true if t is a reference to the nullObject of type T
template<class T>
inline bool isNull(const T& t);

//- Return true if t is not a reference to the nullObject of type T
template<class T>
inline bool notNull(const T& t);


//- Return true if t is a pointer to the nullObject of type T
template<class T>
inline bool isNull(const T* t);

//- Return true if t is not a pointer to the nullObject of type T
template<class T>
inline bool notNull(const T* t);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "nullObjectI.H"

#endif

