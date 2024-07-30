// AN EDGE IS A LIST OF 2-POINT LABELS. THIS SUPPORTS THE DISCRETIZATION 
// ON A 2-D FLAT MESH. CONTAINS VARIOUS MEMBER FUNCTIONS

#ifndef edge_H
#define edge_H

//#include "FixedList.H"
#include "pointField.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class edge;
inline bool operator==(const edge& a, const edge& b);
inline bool operator!=(const edge& a, const edge& b);


/*---------------------------------------------------------------------------*\
                            Class edge Declaration
\*---------------------------------------------------------------------------*/

class edge
:
    public FixedList<label, 2>
{

public:

    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Null constructor for lists
        inline edge();

        //- Construct from components
        inline edge(const label a, const label b);

        //- Construct from FixedList
        inline edge(const FixedList<label, 2>&);

        //- Construct from Istream
        inline edge(Istream&);


    // Member Functions

        //- Return start vertex label
        inline label start() const;

        //- Return start vertex label
        inline label& start();

        //- Return end vertex label
        inline label end() const;

        //- Return end vertex label
        inline label& end();

        //- Return true if connected to given edge
        inline bool connected(const edge& a) const;

        //- Return common vertex
        inline label commonVertex(const edge& a) const;

        //- Given one vertex, return the other
        inline label otherVertex(const label a) const;

        //- Flip the edge in-place.
        inline void flip();

        //- Return reverse edge
        inline edge reverseEdge() const;

        //- Return centre (centroid)
        inline point centre(const pointField&) const;

        //- Return the vector (end - start)
        inline vector vec(const pointField&) const;

        //- Return scalar magnitude
        inline scalar mag(const pointField&) const;

        //- Return edge line
        inline linePointRef line(const pointField&) const;

        //- Compare edges
        //  Returns:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        static inline int compare(const edge&, const edge&);


    // Friend Operators

        friend bool operator==(const edge& a, const edge& b);
        friend bool operator!=(const edge& a, const edge& b);
};


//- Hash specialisation for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e, unsigned seed) const
{
    unsigned val = seed;

    if (e[0] < e[1])
    {
        val = Hash<label>()(e[0], val);
        val = Hash<label>()(e[1], val);
    }
    else
    {
        val = Hash<label>()(e[1], val);
        val = Hash<label>()(e[0], val);
    }

    return val;
}


//- Hash specialisation for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e) const
{
    return Hash<edge>()(e, 0);
}


template<>
inline bool contiguous<edge>()  {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "edgeI.H"

#endif
