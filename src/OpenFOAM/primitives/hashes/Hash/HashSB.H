// HASH FUNCTION CLASS FOR PRIMITIVES.

#ifndef HashSB_H
#define HashSB_H

#include "label.H"
//#include "uLabel.H"
//#include "HasherSB.H"
//#include "pTraits.H"
#include "fileName.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Hash Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class Hash
{
public:

    Hash()
    {}

    unsigned operator()(const PrimitiveType& p, unsigned seed) const
    {
        return Hasher(&p, sizeof(p), seed);
    }

    unsigned operator()(const PrimitiveType& p) const
    {
        return Hasher(&p, sizeof(p));
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hash specialisation for hashing labels
template<>
class Hash<Foam::label>
{
public:

    Hash()
    {}

    //- Incrementally hash a label.
    //  This will necessarily return a different value than the
    //  non-incremental version.
    unsigned operator()(const label p, unsigned seed) const
    {
        return Hasher(&p, sizeof(label), seed);
    }

    //- Return the unsigned representation of a label.
    //  This helps if people have relied on the hash value corresponding to
    //  the natural order.
    unsigned operator()(const label p) const
    {
        return p;
    }
};


//- Hash specialisation for hashing strings
template<>
class Hash<Foam::string>
{
public:

    Hash()
    {}

    unsigned operator()(const string& p, unsigned seed) const
    {
        return string::hash()(p, seed);
    }
    unsigned operator()(const string& p) const
    {
        return string::hash()(p);
    }
};


//- Hash specialisation for hashing words
template<>
class Hash<Foam::word>
{
public:

    Hash()
    {}

    unsigned operator()(const word& p, unsigned seed) const
    {
        return word::hash()(p, seed);
    }
    unsigned operator()(const word& p) const
    {
        return word::hash()(p);
    }
};


//- Hash specialisation for hashing fileNames
template<>
class Hash<Foam::fileName>
{
public:

    Hash()
    {}

    unsigned operator()(const fileName& p, unsigned seed) const
    {
        return fileName::hash()(p, seed);
    }
    unsigned operator()(const fileName& p) const
    {
        return fileName::hash()(p);
    }
};


//- Hash specialisation for hashing wordRes
template<>
class Hash<Foam::wordRe>
{
public:

    Hash()
    {}

    unsigned operator()(const wordRe& p, unsigned seed) const
    {
        return wordRe::hash()(p, seed);
    }
    unsigned operator()(const wordRe& p) const
    {
        return wordRe::hash()(p);
    }
};


//- Hash specialisation for hashing keyTypes
template<>
class Hash<Foam::keyType>
{
public:

    Hash()
    {}

    unsigned operator()(const keyType& p, unsigned seed) const
    {
        return keyType::hash()(p, seed);
    }
    unsigned operator()(const keyType& p) const
    {
        return keyType::hash()(p);
    }
};



//- Hash specialisation for hashing pointer addresses.
//  Treat a pointer like a long.
//  This should work for both 32-bit and 64-bit pointers.
template<>
class Hash<void*>
{
public:

    Hash()
    {}

    unsigned operator()(const void* const& p, unsigned seed) const
    {
        return Hash<long>()(long(p), seed);
    }

    unsigned operator()(const void* const& p) const
    {
        return Hash<long>()(long(p));
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
