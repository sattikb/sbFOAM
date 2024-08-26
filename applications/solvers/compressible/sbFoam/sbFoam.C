// THIS IS THE MAIN SOLVER. ALL BEGIN AND END HERE
//
#include "fvCFDSB.H"
#include "dynamicFvMesh.H"
#include "dynamicMomentumTransportModelSB.H"
#include "fluidThermophysicalTransportModelSB.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    Info<<"SATTIK ENTERING setRootCaseLists.H"<<endl;
    #include "setRootCaseLists.H"
    Info<<"SATTIK EXITED setRootCaseLists.H"<<endl;

    Info<<"SATTIK ENTERING createTime.H"<<endl;
    #include "createTime.H"
    Info<<"SATTIK EXITED createTime.H"<<endl;
    
    Info<<"SATTIK ENTERING createDynamicFvMesh.H"<<endl;
    #include "createDynamicFvMesh.H"
    Info<<"SATTIK EXITED createDynamicFvMesh.H"<<endl;
    
    Info<<"SATTIK ENTERING createDymControls.H"<<endl;
    #include "createDyMControls.H"
    Info<<"SATTIK EXITED createDymControls.H"<<endl;
    
    Info<<"SATTIK ENTERING initContinuityErrors.H"<<endl;
    #include "initContinuityErrs.H"
    Info<<"SATTIK EXITED initContinuityErrors.H"<<endl;
    
    Info<<"SATTIK ENTERING createFields.H"<<endl;
    #include "createFields.H"
    Info<<"SATTIK EXITED createFields.H"<<endl;
    
    Info<<"SATTIK ENTERING createFieldRefs.H"<<endl;
    #include "createFieldRefs.H"
    Info<<"SATTIK EXITED createFieldRefs.H"<<endl;
    
    Info<<"SATTIK ENTERING createRhoUfIfPresent.H"<<endl;
    #include "createRhoUfIfPresent.H"
    Info<<"SATTIK EXITED createRhoUfIfPresent.H"<<endl;

    Info<<"SATTIK VALIDATING TURBULENCE."<<endl;
    turbulence->validate();

    Info<<"SATTIK CHECKING LOCAL TIME STEPPING."<<endl;
    if (!LTS)
    {
    	Info<<"SATTIK HAS LOCAL TIME STEPPING."<<endl;
    	Info<<"SATTIK ENTERING compressibleCourantNo.H"<<endl;
        #include "compressibleCourantNo.H"
    	Info<<"SATTIK EXITED compressibleCourantNo.H"<<endl;
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", rho*U);
                }

                fvModels.preUpdateMesh();

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if
            (
                !mesh.steady()
             && !pimple.simpleRho()
             && pimple.firstPimpleIter()
            )
            {
                #include "rhoEqn.H"
            }

            fvModels.correct();

		//SATTIK VELOCITY DEPENDENT VISCOSITY
		srSB = Foam::sqrt(2.0)*mag(symm(fvc::grad(U)));
		volScalarField srN0SB = max(srSB,nonZeroSmall);
		muSB = muInf + (mu0-muInf)/( 1.0+pow(kSB*srN0SB,nSB) );

            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        if (!mesh.steady())
        {
            rho = thermo.rho();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
