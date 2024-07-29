// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

muOF = thermo.mu();
srSB = Foam::sqrt(2.0)*mag(symm(fvc::grad(U)));
volScalarField srN0SB = max(srSB,nonZeroSmall);
muSB = muOF + muInf + (mu0-muInf)/( 1.0+pow(kSB*srN0SB,nSB) );

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(muSB , U)
  - fvc::div(muSB * dev2(T(fvc::grad(U))))
 ==
    fvModels.source(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
