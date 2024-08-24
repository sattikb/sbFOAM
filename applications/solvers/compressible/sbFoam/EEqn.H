{
    volScalarField& he = thermo.he();

    volScalarField kappa = thermo.kappa();
    volScalarField Cp = thermo.Cp();

    volTensorField StrainRateSB = 0.5*( fvc::grad(U) + T(fvc::grad(U)) );
    volTensorField RotRateSB    = 0.5*( fvc::grad(U) - T(fvc::grad(U)) );

    volTensorField sigma1 = -p*ISB;
    volTensorField sigma2 = -2./3. * muSB * fvc::div(U) * ISB;
    volTensorField sigma3 = muSB*( fvc::grad(U) + T(fvc::grad(U)) );

    volScalarField SBVD1 = fvc::div(sigma1 & U);
    volScalarField SBVD2 = fvc::div(sigma2 & U);
    volScalarField SBVD3 = fvc::div(sigma3 & U);

    SBVD  = SBVD1 + SBVD2 + SBVD3;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt - fvc::div(phi,p/rho)
        )
      - fvm::laplacian(kappa/Cp, he)	//heat conduction term (implicit)
      - ( SBVD1 + SBVD2 + SBVD3 )
     ==
        fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    // CALCULATES THE THERMO/TRANSPORT PROPERTIES USING THE he SOLUTION
    // CAN BE SEEN IN THE heRhoThermo FILE
    thermo.correct();
}
