#include "sbLiqThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

defineTemplateTypeNameAndDebugWithName
(
    sbThermopureMixtureliquidProperties,
    "sbThermo<pureMixture<liquid,sensibleInternalEnergy>>",
    0
);

addToRunTimeSelectionTable
(
    basicThermo,
    sbThermopureMixtureliquidProperties,
    fvMesh
);

addToRunTimeSelectionTable
(
    fluidThermo,
    sbThermopureMixtureliquidProperties,
    fvMesh
);

addToRunTimeSelectionTable
(
    sbRhoThermo,
    sbThermopureMixtureliquidProperties,
    fvMesh
);


defineTemplateTypeNameAndDebugWithName
(
    sbThermopureMixtureEnthalpyliquidProperties,
    "sbThermo<pureMixture<liquid,sensibleEnthalpy>>",
    0
);

addToRunTimeSelectionTable
(
    basicThermo,
    sbThermopureMixtureEnthalpyliquidProperties,
    fvMesh
);

addToRunTimeSelectionTable
(
    fluidThermo,
    sbThermopureMixtureEnthalpyliquidProperties,
    fvMesh
);

addToRunTimeSelectionTable
(
    sbRhoThermo,
    sbThermopureMixtureEnthalpyliquidProperties,
    fvMesh
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
