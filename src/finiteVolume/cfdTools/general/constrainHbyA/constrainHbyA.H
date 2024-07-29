#ifndef constrainHbyA_H
#define constrainHbyA_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<volVectorField> constrainHbyA
(
    const tmp<volVectorField>& tHbyA,
    const volVectorField& U,
    const volScalarField& p
);

tmp<surfaceScalarField> constrainPhiHbyA
(
    const tmp<surfaceScalarField>& tphiHbyA,
    const volVectorField& U,
    const volScalarField& p
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
