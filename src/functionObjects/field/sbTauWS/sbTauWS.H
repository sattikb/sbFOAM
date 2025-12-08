//     Example of function object specification:
//     \verbatim
//     sbTauWS
//     {
//         type        sbTauWS;
//         libs        ("libfieldFunctionObjects.so");
//         ...
//         patches     (".*Wall");
//     }
//     \endverbatim

#ifndef functionObjects_sbTauWS_H
#define functionObjects_sbTauWS_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sbTauWS Declaration
\*---------------------------------------------------------------------------*/

class sbTauWS
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

public:

    //- Runtime type information
    TypeName("sbTauWS");


    // Constructors

        //- Construct from Time and dictionary
        sbTauWS
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        sbTauWS(const sbTauWS&) = delete;


    //- Destructor
    virtual ~sbTauWS();


    // Member Functions

        //- Read the sbTauWS data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sbTauWS&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
