// DIRECTION IS 8-BIT UNSIGNED INT TO REPRESENT CARTESIAN DIRECTIONS

#ifndef directionSB_H
#define directionSB_H

#include <cstdint>
#include <iostream>
#include "charSB.H" 	// TO GET Istream and Ostream

namespace Foam
{

typedef uint8_t direction;

direction readDirection(Istream&);
Istream& operator>>(Istream&, direction&);
Ostream& operator<<(Ostream&, const direction);
std::ostream& operator<<(std::ostream&, const direction);

} // End namespace Foam

#endif
