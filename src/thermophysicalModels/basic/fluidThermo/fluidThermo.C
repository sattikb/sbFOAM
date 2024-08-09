#include "fluidThermoSB.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidThermo, 0);
    defineRunTimeSelectionTable(fluidThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermo::implementation::implementation
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    p_(lookupOrConstruct(mesh, "p"))
{}


Foam::fluidThermo::implementation::implementation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    p_(lookupOrConstruct(mesh, "p"))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidThermo> Foam::fluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<fluidThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermo::~fluidThermo()
{}


Foam::fluidThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}


Foam::volScalarField& Foam::fluidThermo::implementation::p()
{
    return p_;
}


const Foam::volScalarField& Foam::fluidThermo::implementation::p() const
{
    return p_;
}
