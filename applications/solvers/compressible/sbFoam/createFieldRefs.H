const volScalarField& psi = thermo.psi();

//
volScalarField srSB
(
    IOobject
    (
        "srSB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("srSB", dimless/dimTime, 0)
);

volScalarField muSB
(
    IOobject
    (
        "muSB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);


volScalarField SBVD
(
    IOobject
    (
        "VISC_DISS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("VISC DISS",dimensionSet(1,-1,-3,0,0,0,0),1e-20)
);

IOdictionary thermoprops
(
 IOobject
 (
  "thermophysicalProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
 )
);

dimensionedScalar mu0  (thermoprops.lookup("mu0")  );
dimensionedScalar muInf(thermoprops.lookup("muInf"));
dimensionedScalar kSB  (thermoprops.lookup("kSB")  );
dimensionedScalar nSB  (thermoprops.lookup("nSB")  );
dimensionedScalar nonZeroSmall  ("nonZeroSmall",dimensionSet(0,0,-1,0,0,0,0),1.0e-20);
