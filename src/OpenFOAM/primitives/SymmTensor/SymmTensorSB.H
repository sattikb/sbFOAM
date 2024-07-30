//  Templated 3D symmetric tensor derived from VectorSpace adding construction
//  from 6 components, element access using xx(), xy() etc. member functions
//  and the inner-product (dot-product) and outer-product of two Vectors
//  (tensor-product) operators.

#ifndef SymmTensorSB_H
#define SymmTensorSB_H

#include "SphericalTensorSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SymmTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SymmTensor
:
    public VectorSpace<SymmTensor<Cmpt>, Cmpt, 6>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef SymmTensor<label> labelType;


    // Member constants

        //- Rank of SymmTensor is 2
        static const direction rank = 2;


    // Static Data Members

        static const SymmTensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YY, YZ, ZZ };


    // Constructors

        //- Construct null
        inline SymmTensor();

        //- Construct initialised to zero
        inline SymmTensor(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline SymmTensor(const VectorSpace<SymmTensor<Cmpt2>, Cmpt2, 6>&);

        //- Construct given SphericalTensor
        inline SymmTensor(const SphericalTensor<Cmpt>&);

        //- Construct given the six components
        inline SymmTensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
                            const Cmpt tyy, const Cmpt tyz,
                                            const Cmpt tzz
        );

        //- Construct from Istream
        SymmTensor(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& xz() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& yz() const;
            inline const Cmpt& zz() const;

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& xz();
            inline Cmpt& yy();
            inline Cmpt& yz();
            inline Cmpt& zz();

        //- Transpose
        inline const SymmTensor<Cmpt>& T() const;


    // Member Operators

        //- Inherit VectorSpace assignment operators
        using SymmTensor::vsType::operator=;

        //- Assign to given SphericalTensor
        inline void operator=(const SphericalTensor<Cmpt>&);
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 2>
{
public:

    typedef SymmTensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SymmTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
