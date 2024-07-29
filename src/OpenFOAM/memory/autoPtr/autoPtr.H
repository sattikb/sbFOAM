//  An auto-pointer similar to the STL auto_ptr but with automatic casting
//  to a reference to the type and with pointer allocation checking on access.

#ifndef autoPtr_H
#define autoPtr_H

#include <cstddef>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoPtr Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class autoPtr
{
    // Public data

        //- Pointer to object
        mutable T* ptr_;


public:

    typedef T Type;


    // Constructors

        //- Store object pointer
        inline explicit autoPtr(T* = nullptr);

        //- Construct as copy by transferring pointer to this autoPtr and
        //  setting the arguments pointer to nullptr
        inline autoPtr(const autoPtr<T>&);

        //- Construct either by transferring pointer or cloning.
        //  Should only be called with type that supports cloning
        inline autoPtr(const autoPtr<T>&, const bool reuse);


    //- Destructor, delete object if pointer is not nullptr
    inline ~autoPtr();


    // Member Functions

        // Check

            //- Return true if the autoPtr is empty (ie, no pointer set)
            inline bool empty() const;

            //- Return true if the autoPtr valid (ie, the pointer is set)
            inline bool valid() const;


        // Edit

            //- Return object pointer for reuse
            inline T* ptr();

            //- Set pointer to that given.
            //  If object pointer already set issue a FatalError
            inline void set(T*);

            //- If object pointer already set, delete object and set to given
            //  pointer
            inline void reset(T* = nullptr);

            //- Delete object (if the pointer is valid) and set pointer to
            //  nullptr
            inline void clear();


        // Member Operators

            //- Return reference to the object data
            inline T& operator()();

            //- Return const reference to the object data
            inline const T& operator()() const;

            //- Return reference to the object data
            inline T& operator*();

            //- Return const reference to the object data
            inline const T& operator*() const;

            //- Const cast to the underlying type reference
            inline operator const T&() const;

            //- Return object pointer
            inline T* operator->();

            //- Return const object pointer
            inline const T* operator->() const;

            //- Take over the object pointer from parameter
            inline void operator=(T*);

            //- Take over the object pointer from parameter
            inline void operator=(const autoPtr<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "autoPtrI.H"

#endif
