// THE PRIMITIVE FOR LINE AND BASIC CALCULATIONS OF A LINE

#ifndef lineSB_H
#define lineSB_H

#include "PointHitSB.H"
#include "point2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class Istream;
class Ostream;


// Forward declaration of friend functions and operators

template<class Point, class PointRef> class line;

template<class Point, class PointRef>
inline Istream& operator>>(Istream&, line<Point, PointRef>&);

template<class Point, class PointRef>
inline Ostream& operator<<(Ostream&, const line<Point, PointRef>&);


/*---------------------------------------------------------------------------*\
                            Class line Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class line
{
    // Private Data

        PointRef a_, b_;


public:

    // Constructors

        //- Construct from two points
        inline line(const Point& start, const Point& end);

        //- Construct from two points in the list of points
        //  The indices could be from edge etc.
        inline line
        (
            const UList<Point>&,
            const FixedList<label, 2>& indices
        );

        //- Construct from Istream
        inline line(Istream&);


    // Member Functions

        // Access

            //- Return first vertex
            inline PointRef start() const;

            //- Return second vertex
            inline PointRef end() const;


        // Properties

            //- Return centre (centroid)
            inline Point centre() const;

            //- Return scalar magnitude
            inline scalar mag() const;

            //- Return start-end vector
            inline Point vec() const;

            //- Return nearest distance to line from a given point
            //  If the nearest point is on the line, return a hit
            PointHit<Point> nearestDist(const Point& p) const;

            //- Return nearest distance from line to line. Returns distance
            //  and sets both points (one on *this, one on the provided
            //  linePointRef.
            scalar nearestDist
            (
                const line<Point, const Point&>& edge,
                Point& thisPoint,
                Point& edgePoint
            ) const;


    // Ostream operator

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            line&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const line&
        );
};


//- 2D specialisation
template<>
scalar line<point2D, const point2D&>::nearestDist
(
    const line<point2D, const point2D&>& edge,
    point2D& thisPoint,
    point2D& edgePoint
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "lineI.H"

#endif
