// IDENTIFIES PATCH BY NAME, PATCH INDEX, PHYSICAL TYPE, AND INGROUP

#ifndef patchIdentifierSB_H
#define patchIdentifierSB_H

#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

class patchIdentifier;
Ostream& operator<<(Ostream&, const patchIdentifier&);


/*---------------------------------------------------------------------------*\
                       Class patchIdentifier Declaration
\*---------------------------------------------------------------------------*/

class patchIdentifier
{
protected:

    // Protected Data

        //- Name of patch
        word name_;

        //- Index of patch in boundary
        label index_;

        //- Optional physical type
        word physicalType_;

        //- Optional groups patch belongs to
        wordList inGroups_;


public:

    // Constructors

        //- Construct from components
        patchIdentifier
        (
            const word& name,
            const label index,
            const word& physicalType = word::null,
            const wordList& inGroups = wordList()
        );

        //- Construct from dictionary
        patchIdentifier
        (
            const word& name,
            const dictionary&,
            const label index
        );

        //- Construct from geometric patch, resetting the index
        patchIdentifier
        (
            const patchIdentifier&,
            const label index
        );


    //- Destructor
    virtual ~patchIdentifier();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the index of this patch in the boundaryMesh
        label index() const
        {
            return index_;
        }

        //- Return the optional physical type of the patch
        const word& physicalType() const
        {
            return physicalType_;
        }

        //- Return the optional physical type of the patch for modification
        word& physicalType()
        {
            return physicalType_;
        }

        //- Return the optional groups patch belongs to
        const wordList& inGroups() const
        {
            return inGroups_;
        }

        //- Return the optional groups patch belongs to for modification
        wordList& inGroups()
        {
            return inGroups_;
        }

        //- Test if in group
        bool inGroup(const word&) const;

        //- Write patchIdentifier as a dictionary
        void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const patchIdentifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
