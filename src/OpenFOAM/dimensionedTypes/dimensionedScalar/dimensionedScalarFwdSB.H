// FORWARD DECLARES dimensionedScalar WHICH IS A CLASS OF RETURN TYPE SCALAR

#ifndef dimensionedScalarFwdSB_H
#define dimensionedScalarFwdSB_H

#include "scalar.H"

namespace Foam
{

template<class Type>
class dimensioned;

typedef dimensioned<scalar> dimensionedScalar;

} // End namespace Foam

#endif
