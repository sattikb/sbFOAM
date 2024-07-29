#ifndef HasherSB_H
#define HasherSB_H

#include <cstddef>

namespace Foam
{

    //- Bob Jenkins's 96-bit mixer hashing function (lookup3)
    //  \param[in] data - a character stream
    //  \param[in] len  - the number of bytes
    //  \param[in] seed - the previous hash, or an arbitrary value
    unsigned Hasher(const void* data, size_t len, unsigned seed = 0);

} // End namespace Foam

#endif
