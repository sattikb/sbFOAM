// THERMO LOOP HAPPENS HERE. THIS IS WHERE WHEN thermo.correct() WITH
// heRhoThermo IS CALLED IN THE EEqn FILE

#include "heRhoThermo.H"
#include "typeInfo.H"  // For runtime type() checks
#include "word.H"  // For word::null (if not already included via thermo headers)
#include <type_traits>  // For std::enable_if, std::declval
#include "fvc.H"
			
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Trait to check if mu3ArgSB exists on the type
template<class T>
struct has_mu3ArgSB {
private:
    template<class U>
    static auto test(int) -> decltype(std::declval<U>().mu3ArgSB(std::declval<Foam::scalar>(), std::declval<Foam::scalar>(), std::declval<Foam::scalar>()), std::true_type{});
    template<class> static std::false_type test(...);
public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

// SFINAE helper: Use mu3ArgSB if available, else standard mu(p, T)
namespace Foam
{
    // Primary overload: Enabled only if has_mu3ArgSB<T>::value == true
    template<class TransportType>
    auto callMu3ArgSB(const TransportType& t, scalar p, scalar T, scalar sr)
        -> typename std::enable_if<has_mu3ArgSB<TransportType>::value, scalar>::type
    {
        return t.mu3ArgSB(p, T, sr);
    }

    // Fallback overload: Enabled only if has_mu3ArgSB<T>::value == false
    template<class TransportType>
    auto callMu3ArgSB(const TransportType& t, scalar p, scalar T, scalar sr)
        -> typename std::enable_if<!has_mu3ArgSB<TransportType>::value, scalar>::type
    {
        return t.mu(p, T);  // Ignores umag
    }
}


template<class BasicRhoThermo, class MixtureType>
void Foam::heRhoThermo<BasicRhoThermo, MixtureType>::calculate()
{
    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    
    //
    bool useStrainRate = this->db().template foundObject<volVectorField>("U");
    dimensionedScalar smallSrSB("smallSrSB", dimless/dimTime, 1e-15);
    volScalarField srN0SB
    (
        IOobject
        (
            "srN0SB",
            this->db().time().timeName(),
            this->db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->p().mesh(),
	smallSrSB
    );

    if (useStrainRate)
    {
        const volVectorField& U = this->db().template lookupObject<volVectorField>("U");
        volTensorField srTensor = (fvc::grad(U) + T(fvc::grad(U)));
        volScalarField srSB = Foam::sqrt(srTensor && fvc::grad(U));
        srN0SB = max(srSB, smallSrSB);
    }

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

        CpCells[celli] = thermoMixture.Cp(pCells[celli], TCells[celli]);
        CvCells[celli] = thermoMixture.Cv(pCells[celli], TCells[celli]);
        psiCells[celli] = thermoMixture.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = thermoMixture.rho(pCells[celli], TCells[celli]);

	/////////////
	muCells[celli] = Foam::callMu3ArgSB(transportMixture, pCells[celli], TCells[celli], srN0SB[celli]);
	/////////////
        // muCells[celli] = transportMixture.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] =
            transportMixture.kappa(pCells[celli], TCells[celli])
           /thermoMixture.Cp(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& CpBf =
        this->Cp_.boundaryFieldRef();

    volScalarField::Boundary& CvBf =
        this->Cv_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

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

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);
		////////////// NEW: Conditional mu call (mirrors internal)
		pmu[facei] = Foam::callMu3ArgSB(transportMixture, pp[facei], pT[facei], srN0SB.boundaryField()[patchi][facei]);
                //////////////
        //        pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cp(pp[facei], pT[facei]);
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

                pCp[facei] = thermoMixture.Cp(pp[facei], pT[facei]);
                pCv[facei] = thermoMixture.Cv(pp[facei], pT[facei]);
                ppsi[facei] = thermoMixture.psi(pp[facei], pT[facei]);
                prho[facei] = thermoMixture.rho(pp[facei], pT[facei]);

		////////////// NEW: Conditional mu call (mirrors internal)
		pmu[facei] = Foam::callMu3ArgSB(transportMixture, pp[facei], pT[facei], srN0SB.boundaryField()[patchi][facei]);
                //////////////
        //        pmu[facei] = transportMixture.mu(pp[facei], pT[facei]);
                palpha[facei] =
                    transportMixture.kappa(pp[facei], pT[facei])
                   /thermoMixture.Cp(pp[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicRhoThermo, class MixtureType>
Foam::heRhoThermo<BasicRhoThermo, MixtureType>::heRhoThermo
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
Foam::heRhoThermo<BasicRhoThermo, MixtureType>::~heRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// THIS IS USUALLY CALLED IN EEqn.H FILE VIA thermo.correct()
template<class BasicRhoThermo, class MixtureType>
void Foam::heRhoThermo<BasicRhoThermo, MixtureType>::correct()
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
