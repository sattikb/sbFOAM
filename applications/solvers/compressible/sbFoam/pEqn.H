if ((!mesh.steady() && !pimple.simpleRho()) || pimple.consistent())
{
    Info<<"SATTIK (UNSTEADY & NOT SIMPLE) OR INCONSITSTENT"<<endl;
    Info<<"	THUS WILL USE EOS TO CALCULATE RHO NOW"<<endl;
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

// THINK OF AU-H=-grad(p)
// THIS IS 1/A -> A is the diagonal matrix from MU=-grad(p)
const volScalarField rAU("rAU", 1.0/UEqn.A());
// THIS IS rho*inv(A)
const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
// THIS IS 1(1/A-H1)
tmp<volScalarField> rAtU
(
    pimple.consistent()
  ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
  : tmp<volScalarField>(nullptr)
);
// THIS IS rho*inv(A-H1)
tmp<surfaceScalarField> rhorAtUf
(
    pimple.consistent()
  ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(rho*rAtU()))
  : tmp<surfaceScalarField>(nullptr)
);
// THIS IS EITHER inv(A-H1) OR inv(A)
const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
// THIS IS EITHER rho*inv(A-H1) OR rho*inv(A)
const surfaceScalarField& rhorAAtUf =
    pimple.consistent() ? rhorAtUf() : rhorAUf;
// THIS is inv(A)&H
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPiso() <= 1)
{
    Info<<"SATTIK  HAS nCorrPiso LESS THAN 1"<<endl;
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
  + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
);

fvc::makeRelative(phiHbyA, rho, U);
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

bool adjustMass = false;

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAAtUf, MRF);

if (pimple.transonic())
{
    Info<<"SATTIK HAS TRANSONIC FLOW"<<endl;
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/fvc::interpolate(rho);

    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(p);
    }

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA) + fvm::div(phid, p)
     ==
        fvModels.source(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    Info<<"SATTIK DOES NOT HAVE TRANSONIC FLOW"<<endl;
    if (mesh.steady())
    {
	Info<<"SATTIK HAS STEADY MESH"<<endl;
        adjustMass = adjustPhi(phiHbyA, U, p);
    }

    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(p)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(p);
    }

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p))
      + fvc::div(phiHbyA)
     ==
        fvModels.source(psi, p, rho.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, p));

        pEqn.setReference
        (
            pressureReference.refCell(),
            pressureReference.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

if (mesh.steady())
{
    #include "incompressible/continuityErrs.H"
}
else
{
    const bool constrained = fvConstraints.constrain(p);

    // Thermodynamic density update
    thermo.correctRho(psi*p - psip0);

    if (constrained)
    {
        rho = thermo.rho();
    }

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvConstraints.constrain(U);
K = 0.5*magSqr(U);

if (mesh.steady())
{
    fvConstraints.constrain(p);
}

// For steady compressible closed-volume cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass && !thermo.incompressible())
{
    p += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(psi);
    p.correctBoundaryConditions();
}

if (mesh.steady() || pimple.simpleRho() || adjustMass)
{
    rho = thermo.rho();
}

// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, rho, U, phi);

if ((mesh.steady() || pimple.simpleRho()) && !pimple.transonic())
{
    rho.relax();
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
