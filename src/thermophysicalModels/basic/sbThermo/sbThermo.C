// THERMO LOOP HAPPENS HERE. THIS IS WHERE WHEN thermo.correct() WITH
// sbThermo IS CALLED IN THE EEqn FILE

#include "sbThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
void Foam::sbThermo<BasicRhoThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();

    Info<<"SATTIK REACHED HERE"<<endl;

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoMixtureType& thermoMixture =
            this->cellThermoMixture(celli);

        const typename MixtureType::transportMixtureType& transportMixture =
            this->cellTransportMixture(celli, thermoMixture);

	// .THE() IS DEFINED IN thermoI.H	
	// CALCULATES T FOR NEW TIME-STEP
        TCells[celli] = thermoMixture.THE
        (
            hCells[celli],	//he SOLUTION FROM SOLVER
            pCells[celli],	// P SOLUTION FROM SOLVER
            TCells[celli]	// T FROM PEVIOUS TIME STEP
        );

    }

    volScalarField::Boundary& pBf = this->p_.boundaryFieldRef();
    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& heBf = this->he().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                phe[facei] = thermoMixture.HE(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoMixtureType& thermoMixture =
                    this->patchFaceThermoMixture(patchi, facei);

                const typename MixtureType::transportMixtureType&
                    transportMixture =
                    this->patchFaceTransportMixture
                    (patchi, facei, thermoMixture);

                pT[facei] = thermoMixture.THE(phe[facei], pp[facei], pT[facei]);

            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::sbThermo<BasicRhoThermo, MixtureType>::sbThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicRhoThermo, MixtureType>(mesh, phaseName)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::sbThermo<BasicRhoThermo, MixtureType>::~sbThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// THIS IS USUALLY CALLED IN EEqn.H FILE VIA thermo.correct()
template<class BasicRhoThermo, class MixtureType>
void Foam::sbThermo<BasicRhoThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();	//CALCULATE IS DEFINED ABOVE

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}
