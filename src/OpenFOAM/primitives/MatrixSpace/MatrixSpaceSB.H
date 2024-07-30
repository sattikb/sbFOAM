// TEMPLATE OF MATRIX SPACE WHOSE ARGUMENTS ARE THE Form OF THE MATRIX 
// SPACE, WHICH WILL FURTHER BE USED TO CREATE THE TUPE OF ELEMENTS, AND NO 
// OF ROWS AND COLUMNS

#ifndef MatrixSpaceSB_H
#define MatrixSpaceSB_H

#include "VectorSpaceSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MatrixSpace Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Cmpt, direction Mrows, direction Ncols>
class MatrixSpace
:
    public VectorSpace<Form, Cmpt, Mrows*Ncols>
{

public:

    //- MatrixSpace type
    typedef MatrixSpace<Form, Cmpt, Mrows, Ncols> msType;


    // Member constants

        static const direction mRows = Mrows;
        static const direction nCols = Ncols;


    // Static member functions

        //- Return the number of rows
        static direction m()
        {
            return Mrows;
        }

        //- Return the number of columns
        static direction n()
        {
            return Ncols;
        }

        //- Return the identity matrix for square matrix spaces
        inline static msType identity();


    // Sub-Block Classes

        //- Const sub-block type
        template<class SubTensor, direction BRowStart, direction BColStart>
        class ConstBlock
        {
            //- Reference to parent matrix
            const msType& matrix_;

        public:

            static const direction mRows = SubTensor::mRows;
            static const direction nCols = SubTensor::nCols;

            //- Return the number of rows in the block
            static direction m()
            {
                return mRows;
            }

            //- Return the number of columns in the block
            static direction n()
            {
                return nCols;
            }

            //- Construct for the given matrix
            inline ConstBlock(const msType& matrix);

            //- Construct and return the sub-ensor corresponding to this block
            inline SubTensor operator()() const;

            //- (i, j) const element access operator
            inline const Cmpt& operator()
            (
                const direction i,
                const direction j
            ) const;
        };


        //- Sub-block type
        template
        <
            class SubTensor,
            direction BRowStart,
            direction BColStart
        >
        class Block
        {
            //- Reference to parent matrix
            msType& matrix_;

        public:

            static const direction mRows = SubTensor::mRows;
            static const direction nCols = SubTensor::nCols;

            //- Return the number of rows in the block
            static direction m()
            {
                return mRows;
            }

            //- Return the number of columns in the block
            static direction n()
            {
                return nCols;
            }

            //- Construct for the given matrix
            inline Block(msType& matrix);

            //- Assignment to a matrix
            template<class Form2>
            inline void operator=
            (
                const MatrixSpace
                <
                    Form2,
                    Cmpt,
                    SubTensor::mRows,
                    SubTensor::nCols
                >& matrix
            );

            //- Assignment to a column vector
            template<class VSForm>
            inline void operator=
            (
                const VectorSpace<VSForm, Cmpt, SubTensor::mRows>& v
            );

            //- Construct and return the sub-tensor corresponding to this block
            inline SubTensor operator()() const;

            //- (i, j) const element access operator
            inline const Cmpt& operator()
            (
                const direction i,
                const direction j
            ) const;

            //- (i, j) element access operator
            inline Cmpt& operator()(const direction i, const direction j);
        };


    // Constructors

        //- Construct null
        inline MatrixSpace();

        //- Construct initialised to zero
        inline MatrixSpace(const Foam::zero);

        //- Construct as copy of a VectorSpace with the same size
        template<class Form2, class Cmpt2>
        inline explicit MatrixSpace
        (
            const VectorSpace<Form2, Cmpt2, Mrows*Ncols>&
        );

        //- Construct from a block of another matrix space
        template
        <
            template<class, direction, direction> class Block2,
            direction BRowStart,
            direction BColStart
        >
        inline MatrixSpace
        (
            const Block2<Form, BRowStart, BColStart>& block
        );

        //- Construct from Istream
        MatrixSpace(Istream&);


    // Member Functions

        //- Fast const element access using compile-time addressing
        template<direction Row, direction Col>
        inline const Cmpt& elmt() const;

        //- Fast element access using compile-time addressing
        template<direction Row, direction Col>
        inline Cmpt& elmt();

        // Const element access functions for a 3x3
        // Compile-time errors are generated for inappropriate use

            inline const Cmpt& xx() const;
            inline const Cmpt& xy() const;
            inline const Cmpt& xz() const;
            inline const Cmpt& yx() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& yz() const;
            inline const Cmpt& zx() const;
            inline const Cmpt& zy() const;
            inline const Cmpt& zz() const;

        // Element access functions for a 3x3
        // Compile-time errors are generated for inappropriate use

            inline Cmpt& xx();
            inline Cmpt& xy();
            inline Cmpt& xz();
            inline Cmpt& yx();
            inline Cmpt& yy();
            inline Cmpt& yz();
            inline Cmpt& zx();
            inline Cmpt& zy();
            inline Cmpt& zz();

        //- Return the transpose of the matrix
        inline typename typeOfTranspose<Cmpt, Form>::type T() const;

        //- Return a const sub-block corresponding to the specified type
        //  starting at the specified row and column
        template<class SubTensor, direction BRowStart, direction BColStart>
        inline ConstBlock<SubTensor, BRowStart, BColStart> block() const;

        //- Return a sub-block corresponding to the specified type
        //  starting at the specified row and column
        template<class SubTensor, direction BRowStart, direction BColStart>
        inline Block<SubTensor, BRowStart, BColStart> block();

        //- (i, j) const element access operator
        inline const Cmpt& operator()
        (
            const direction& i,
            const direction& j
        ) const;

        //- (i, j) element access operator
        inline Cmpt& operator()(const direction& i, const direction& j);


    // Member Operators

        //- Assignment to zero
        inline void operator=(const Foam::zero);

        //- Assignment to a block of another matrix space
        template
        <
            template<class, direction, direction> class Block2,
            direction BRowStart,
            direction BColStart
        >
        inline void operator=
        (
            const Block2<Form, BRowStart, BColStart>& block
        );

        //- Inner product with a compatible square matrix
        template<class Form2>
        inline void operator&=
        (
            const MatrixSpace<Form, Cmpt, Ncols, Ncols>& matrix
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "MatrixSpaceI.H"

#endif
