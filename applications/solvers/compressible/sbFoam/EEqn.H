{
    volScalarField& he = thermo.he();

    volScalarField kappa = thermo.kappa();
    volScalarField Cp = thermo.Cp();
	srSB = Foam::sqrt(2.0)*mag(symm(fvc::grad(U)));
	volScalarField srN0SB = max(srSB,nonZeroSmall);
	muSB = muInf + (mu0-muInf)/( 1.0+pow(kSB*srN0SB,nSB) );

    volTensorField TauSB = muSB * dev2(T(fvc::grad(U)));

	volTensorField sigmaSB = -p*ISB + muSB*symm(fvc::grad(U) );

//    SBVD = TauSB && fvc::grad(U);
//    SBVD = fvc::div(TauSB & U);
    SBVD = fvc::div(sigmaSB & U);

//    volScalarField sattikViscDiss = fvc::div(muSB * (TauSB & U));

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
      - SBVD 
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
