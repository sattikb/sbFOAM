#include "dynamicMomentumTransportModelsSB.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeBaseMomentumTransportModel
(
    geometricOneField,
    volScalarField,
    compressibleMomentumTransportModel,
    CompressibleMomentumTransportModel,
    dynamicTransportModel
);


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "generalisedNewtonian.H"
makeLaminarModel(generalisedNewtonian);

#include "lambdaThixotropic.H"
makeLaminarModel(lambdaThixotropic);

#include "Maxwell.H"
makeLaminarModel(Maxwell);

#include "Giesekus.H"
makeLaminarModel(Giesekus);

#include "PTT.H"
makeLaminarModel(PTT);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "SpalartAllmaras.H"
makeRASModel(SpalartAllmaras);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "RNGkEpsilon.H"
makeRASModel(RNGkEpsilon);

#include "realizableKE.H"
makeRASModel(realizableKE);

#include "buoyantKEpsilon.H"
makeRASModel(buoyantKEpsilon);

#include "LaunderSharmaKE.H"
makeRASModel(LaunderSharmaKE);

#include "kOmega.H"
makeRASModel(kOmega);

#include "kOmega2006.H"
makeRASModel(kOmega2006);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "kOmegaSSTSAS.H"
makeRASModel(kOmegaSSTSAS);

#include "kOmegaSSTLM.H"
makeRASModel(kOmegaSSTLM);

#include "v2f.H"
makeRASModel(v2f);

#include "LRR.H"
makeRASModel(LRR);

#include "SSG.H"
makeRASModel(SSG);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "WALE.H"
makeLESModel(WALE);

#include "kEqn.H"
makeLESModel(kEqn);

#include "dynamicKEqn.H"
makeLESModel(dynamicKEqn);

#include "dynamicLagrangian.H"
makeLESModel(dynamicLagrangian);

#include "kOmegaSSTDES.H"
makeLESModel(kOmegaSSTDES);

#include "SpalartAllmarasDES.H"
makeLESModel(SpalartAllmarasDES);

#include "SpalartAllmarasDDES.H"
makeLESModel(SpalartAllmarasDDES);

#include "SpalartAllmarasIDDES.H"
makeLESModel(SpalartAllmarasIDDES);

#include "DeardorffDiffStress.H"
makeLESModel(DeardorffDiffStress);


// ************************************************************************* //
