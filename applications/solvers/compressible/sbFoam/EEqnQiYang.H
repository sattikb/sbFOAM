{
    volScalarField& he = thermo.he();

    volScalarField kappa = thermo.kappa();
    volScalarField Cp = thermo.Cp();

    // FROM reactingFoam Yang //////////////////////
    volTensorField sigma1 = -2./3. * muSB * fvc::div(U) * ISB;
    volTensorField sigma2 = muSB*( fvc::grad(U) + T(fvc::grad(U)) );

    volScalarField SBVD1 = sigma1 && fvc::grad(U);
    volScalarField SBVD2 = sigma2 && fvc::grad(U);

    volScalarField UdotgradP = U & fvc::grad(p);
    ////////////////////////////////////////////

    SBVD  = SBVD2 + SBVD2;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt - UdotgradP
        )
     ==
        fvm::laplacian(kappa/Cp, he)	//heat conduction term (implicit)
      + SBVD
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    // CALCULATES THE THERMO/TRANSPORT PROPERTIES USING THE he SOLUTION
    // CAN BE SEEN IN THE heRhoThermo FILE
    thermo.correct();
}
