// CALCULATES AND OUTPUTS THE MEAN AND MAX Co Nos
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()/rho.primitiveField()
    );

    Info<<"SATTIK VALUES: "<<runTime.deltaTValue()<<endl;

    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum << " max: " << CoNum << endl;

