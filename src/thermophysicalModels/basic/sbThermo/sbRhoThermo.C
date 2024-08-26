#include "sbRhoThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sbRhoThermo, 0);
    defineRunTimeSelectionTable(sbRhoThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sbRhoThermo::implementation::implementation
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


Foam::sbRhoThermo::implementation::implementation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    rho_
    (
        IOobject
        (
            phasePropertyName("thermo:rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimDensity
    ),

    psi_
    (
        IOobject
        (
            phasePropertyName("thermo:psi", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, -2, 2, 0, 0)
    ),

    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sbRhoThermo> Foam::sbRhoThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<sbRhoThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sbRhoThermo::~sbRhoThermo()
{}


Foam::sbRhoThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::sbRhoThermo::implementation::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::sbRhoThermo::implementation::rho
(
    const label patchi
) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::sbRhoThermo::implementation::rho()
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::sbRhoThermo::implementation::rho0() const
{
    return rho_.oldTime();
}


void Foam::sbRhoThermo::implementation::correctRho(const volScalarField& deltaRho)
{
    rho_ += deltaRho;
}


const Foam::volScalarField& Foam::sbRhoThermo::implementation::psi() const
{
    return psi_;
}


Foam::tmp<Foam::volScalarField> Foam::sbRhoThermo::implementation::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::sbRhoThermo::implementation::mu
(
    const label patchi
) const
{
    return mu_.boundaryField()[patchi];
}
