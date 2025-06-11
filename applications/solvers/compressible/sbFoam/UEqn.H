// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

// FROM KENNETH KUO ///////////////////////
volTensorField tau1 = -2./3. * muSB * fvc::div(U) * ISB;
volTensorField tau2 = muSB*( fvc::grad(U) + T(fvc::grad(U)) );

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(muSB , U)	//same as div(mu*grad(U))
 ==
   fvc::div(muSB * T(fvc::grad(U)))
  + fvc::div(tau1) 
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    Info<<"SATTIK GOING FOR U FVCONSTRAINTS IN MOMENTUM PREDICTOR LOOP"<<endl;
    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
