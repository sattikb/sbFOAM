//  Templated 3D Vector derived from VectorSpace adding construction from
//  3 components, element access using x(), y() and z() member functions and
//  the inner-product (dot-product) and cross product operators.
//
//  A centre() member function which returns the Vector for which it is called
//  is defined so that point which is a typedef to Vector\<scalar\> behaves as
//  other shapes in the shape hierarchy.

#ifndef VectorSB_H
#define VectorSB_H

#include "VectorSpaceSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T> class List;

/*---------------------------------------------------------------------------*\
                           Class Vector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Vector
:
    public VectorSpace<Vector<Cmpt>, Cmpt, 3>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Vector<label> labelType;


    // Member constants

        //- Rank of Vector is 1
        static const direction rank = 1;


    //- Component labeling enumeration
    enum components { X, Y, Z };


    // Constructors

        //- Construct null
        inline Vector();

        //- Construct initialised to zero
        inline Vector(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline Vector(const VectorSpace<Vector<Cmpt2>, Cmpt2, 3>&);

        //- Construct given three components
        inline Vector(const Cmpt& vx, const Cmpt& vy, const Cmpt& vz);

        //- Construct from Istream
        inline Vector(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& x() const;
            inline const Cmpt& y() const;
            inline const Cmpt& z() const;

            inline Cmpt& x();
            inline Cmpt& y();
            inline Cmpt& z();

        //- Return *this (used for point which is a typedef to Vector<scalar>.
        inline const Vector<Cmpt>& centre
        (
            const Foam::List<Vector<Cmpt>>&
        ) const;
};


template<class Cmpt>
class typeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "VectorI.H"

#endif
