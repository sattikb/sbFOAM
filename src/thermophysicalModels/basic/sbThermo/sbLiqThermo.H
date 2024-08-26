#include "sbRhoThermo.H"
#include "sbThermo.H"
#include "pureMixture.H"
#include "thermo.H"
#include "sensibleInternalEnergy.H"
#include "sensibleEnthalpySB.H"
#include "thermophysicalPropertiesSelector.H"
#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef sbThermo
<
    sbRhoThermo::composite,
    pureMixture
    <
        species::thermo
        <
            thermophysicalPropertiesSelector<liquidProperties>,
            sensibleInternalEnergy
        >
    >
> sbThermopureMixtureliquidProperties;

typedef sbThermo
<
    sbRhoThermo::composite,
    pureMixture
    <
        species::thermo
        <
            thermophysicalPropertiesSelector<liquidProperties>,
            sensibleEnthalpy
        >
    >
> sbThermopureMixtureEnthalpyliquidProperties;

} // End namespace Foam
