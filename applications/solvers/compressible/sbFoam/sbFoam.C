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

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
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
	    volTensorField srTensor = (fvc::grad(U) + T(fvc::grad(U)));
	    srSB = Foam::sqrt(srTensor && fvc::grad(U));
	    volScalarField srN0SB = max(srSB,nonZeroSmall);
	    muSB = muInf + (mu0-muInf)/( 1.0+pow(kSB*srN0SB,nSB) );
	    //NEAR WALL MODIFICATIONS
	    forAll(mesh.boundary(), patchI)
	    {
		const fvPatch& patch = mesh.boundary()[patchI];
		const labelList& faceCells = patch.faceCells(); // Cells attached to patch

		forAll(faceCells, i)
		{
		    label cellID = faceCells[i];

		    if (mesh.cellCells()[cellID].size() > 2) // Ensure three layers exist
		    {
			muWall[cellID] = muSB[cellID] * 2;
			muWall[mesh.cellCells()[cellID][0]] = muSB[mesh.cellCells()[cellID][0]] * 20;
			muWall[mesh.cellCells()[cellID][1]] = muSB[mesh.cellCells()[cellID][1]] * 20;
		    }
		}
	    }
	   // muSB = muSB + muWall;
		

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
		Info<<"SATTIK IN TURBCORRECTOR LOOP."<<endl;
                turbulence->correct();
                thermophysicalTransport->correct();
            }
	    
            #include "EEqn.H"
           // #include "TEqn.H"
            #include "computeTauRT.H"

	//   forAll(tauRT, celli)
	//	{
	//	    const vector& C = mesh.C()[celli]; // Cell center coordinates
	//	    tauRT[celli] = 6654212.9 / (C.x() * C.x() + C.y() * C.y());
	//	}

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
