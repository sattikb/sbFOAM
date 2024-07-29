#include "constrainHbyA.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::constrainHbyA
(
    const tmp<volVectorField>& tHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<volVectorField> tHbyANew;

    if (tHbyA.isTmp())
    {
        tHbyANew = tHbyA;
        tHbyANew.ref().rename(IOobject::groupName("HbyA", U.group()));
    }
    else
    {
        tHbyANew = volVectorField::New
        (
            IOobject::groupName("HbyA", U.group()),
            tHbyA
        );
    }

    volVectorField& HbyA = tHbyANew.ref();
    volVectorField::Boundary& HbyAbf = HbyA.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            HbyAbf[patchi] = U.boundaryField()[patchi];
        }
    }

    return tHbyANew;
}


Foam::tmp<Foam::surfaceScalarField> Foam::constrainPhiHbyA
(
    const tmp<surfaceScalarField>& tphiHbyA,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<surfaceScalarField> tphiHbyANew;

    if (tphiHbyA.isTmp())
    {
        tphiHbyANew = tphiHbyA;
        tphiHbyANew.ref().rename(IOobject::groupName("phiHbyA", U.group()));
    }
    else
    {
        tphiHbyANew = surfaceScalarField::New
        (
            IOobject::groupName("phiHbyA", U.group()),
            tphiHbyA
        );
    }

    surfaceScalarField& phiHbyA = tphiHbyANew.ref();
    surfaceScalarField::Boundary& phiHbyAbf = phiHbyA.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            phiHbyAbf[patchi] =
                U.mesh().Sf().boundaryField()[patchi]
              & U.boundaryField()[patchi];
        }
    }

    return tphiHbyANew;
}


// ************************************************************************* //
