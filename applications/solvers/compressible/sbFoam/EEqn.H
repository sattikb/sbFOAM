{
    volScalarField& he = thermo.he();

    volScalarField muSB = thermo.mu();
    volScalarField kappa = thermo.kappa();
    volScalarField Cp = thermo.Cp();

    ////        SHEAR STRESS TAU
    volTensorField tau1 = -2./3. * muSB * fvc::div(U) * ISB;
    volTensorField tau2 = muSB*( fvc::grad(U) + T(fvc::grad(U)) );

    volScalarField SBVD1 = fvc::div(tau1 & U);
    volScalarField SBVD2 = fvc::div(tau2 & U);

    volScalarField SBVD  = SBVD1 + SBVD2;
    ////

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
     ==
        fvModels.source(rho, he)
      + fvm::laplacian(kappa/Cp, he)
      + SBVD
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();
}
