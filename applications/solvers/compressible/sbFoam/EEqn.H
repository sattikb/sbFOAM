{
    Info<<"SATTIK IN  ENERGY EQUATION FILE"<<endl;
    volScalarField& he = thermo.he();
    heSol = thermo.he();

    volScalarField kappa = thermo.kappa();
    volScalarField Cp = thermo.Cp();

    ////	SHEAR STRESS TAU
    volTensorField tau1 = -2./3. * muSB * fvc::div(U) * ISB;
    volTensorField tau2 = muSB*( fvc::grad(U) + T(fvc::grad(U)) );
    volTensorField tauTot = tau1 + tau2;

    volScalarField SBVD1 = fvc::div(tau1 & U);
    volScalarField SBVD2 = fvc::div(tau2 & U);

    SBVD  = SBVD1 + SBVD2;
    ////
    Info<<"SATTIK IN  ENERGY EQUATION FILE"<<endl;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -fvc::ddt(p)
        )
     ==
        fvm::laplacian(kappa/Cp, he)	//implicit heat conduction term
      + SBVD
    );

    Info<<"SATTIK RELAXING ENERGY EQUATION"<<endl;
    EEqn.relax();

    fvConstraints.constrain(EEqn);

    Info<<"SATTIK SOLVING ENERGY EQUATION"<<endl;
    EEqn.solve();

    fvConstraints.constrain(he);

    // CALCULATES THE THERMO/TRANSPORT PROPERTIES USING THE he SOLUTION
    // CAN BE SEEN IN THE heRhoThermo FILE

    Info<<"SATTIK CORRECTING ENERGY EQUATION"<<endl;
    thermo.correct();

}
