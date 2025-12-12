#include "sbTauWS.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "kinematicMomentumTransportModel.H"
#include "dynamicMomentumTransportModelSB.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sbTauWS, 0);
    addToRunTimeSelectionTable(functionObject, sbTauWS, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::sbTauWS::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall shear stress");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sbTauWS::sbTauWS
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sbTauWS::~sbTauWS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sbTauWS::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    resetName(typeName);
    resetLocalObjectName(typeName);

    return true;
}


bool Foam::functionObjects::sbTauWS::execute()
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const volScalarField& mu = mesh_.lookupObject<volScalarField>("thermo:mu");
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();

    tmp<volVectorField> ttau
    (
        new volVectorField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimensionSet(1, -1, -1, 0, 0), Zero)
        )
    );

    volVectorField::Boundary& tauBf = ttau.ref().boundaryFieldRef();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const label patchi = iter.key();
        const fvPatch& pPatch = mesh_.boundary()[patchi];

        const scalarField& mup = mu.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const tensorField gradUp = gradU.boundaryField()[patchi];

        const vectorField np = pPatch.nf();

	tensorField shearTensor = mup*( gradUp + gradUp.T() - 2.0/3.0*tr(gradUp)*I );
	tensorField cauchyTensor = shearTensor - pp*I;
        
	vectorField traction = cauchyTensor & np;
	vectorField tauWSVal = traction - (traction & np)*np;

//	const vectorField& Uf = U.boundaryField()[patchi];       
//        const vectorField  Uc = U.boundaryField()[patchi].patchInternalField();
//	const scalarField delta = mesh_.nonOrthDeltaCoeffs().boundaryField()[patchi];
//	vectorField tauHat = tauWSVal/mag(tauWSVal);
//	vectorField UcTan = (Uc & tauHat)*tauHat;
//        vectorField tauWSVal2 = mup * (Uf - UcTan) * delta;   

        tauBf[patchi] = tauWSVal;
    }

    // Store the field so it can be written
    return store(type(), ttau);
}


bool Foam::functionObjects::sbTauWS::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volVectorField& sbTauWS =
        obr_.lookupObject<volVectorField>(type());

    const fvPatchList& patches = mesh_.boundary();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = sbTauWS.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            file() << mesh_.time().value()
                << tab << pp.name()
                << tab << minSsp
                << tab << maxSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
