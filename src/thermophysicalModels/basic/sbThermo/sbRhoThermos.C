#include "sbRhoThermo.H"
#include "sbThermo.H"
#include "pureMixture.H"

#include "forGases.H"
#include "forLiquids.H"
#include "forTabulated.H"
#include "makeThermo.H"

namespace Foam
{
    forGases(makeThermo, sbRhoThermo, sbThermo, pureMixture);
    forLiquids(makeThermo, sbRhoThermo, sbThermo, pureMixture);
    forTabulated(makeThermo, sbRhoThermo, sbThermo, pureMixture);
}
