{
    volScalarField kappa = thermo.kappa();
    volScalarField Cp = thermo.Cp();

    ////	SHEAR STRESS TAU
    volTensorField tau1 = -2./3. * muSB * fvc::div(U) * ISB;
    volTensorField tau2 = muSB*( fvc::grad(U) + T(fvc::grad(U)) );
    volTensorField tauTot = tau1 + tau2;

    volScalarField SBVD1 = fvc::div(tau1 & U);
    volScalarField SBVD2 = fvc::div(tau2 & U);

    SBVD  = tauTot && fvc::grad(U);
    ////

    fvScalarMatrix TSBEqn
    (
        Cp*fvm::ddt(rho, TSB)
       	+ Cp*fvm::div(phi, TSB)
     ==
        fvm::laplacian(kappa, TSB)	//implicit heat conduction term
      + SBVD
    );

    TSBEqn.solve();

}
