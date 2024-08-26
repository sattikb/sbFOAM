// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);


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

Info<<"SATTIK GOING FOR U FVCONSTRAINTS IN UEqn.H"<<endl;
fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    Info<<"SATTIK GOING FOR U FVCONSTRAINTS IN MOMENTUM PREDICTOR LOOP"<<endl;
    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
