// CLASS FOR HANDLING WORDS DERIVED FROM STRINGS.

#ifndef word_H
#define word_H

#include "stringSB.H"

namespace Foam
{

// Forward declaration of friend functions and operators
class word;
inline word operator&(const word&, const word&);
Istream& operator>>(Istream&, word&);
Ostream& operator<<(Ostream&, const word&);


/*---------------------------------------------------------------------------*\
                            Class word Declaration
\*---------------------------------------------------------------------------*/

class word
:
    public string
{
    // Private Member Functions

        //- Strip invalid characters from this word
        inline void stripInvalid();


public:

    // Static Data Members

        static const char* const typeName;
        static int debug;

        //- An empty word
        static const word null;


    // Constructors

        //- Construct null
        inline word();

        //- Copy constructor
        inline word(const word&);

        //- Copy constructor of character array
        inline word(const char*, const bool doStripInvalid=true);

        //- Copy constructor with a maximum number of characters
        inline word
        (
            const char*,
            const size_type,
            const bool doStripInvalid
        );

        //- Copy constructor of string
        inline word(const string&, const bool doStripInvalid=true);

        //- Copy constructor of std::string
        inline word(const std::string&, const bool doStripInvalid=true);

        //- Construct from Istream
        word(Istream&);


    // Member Functions

        //- Is this character valid for a word
        inline static bool valid(char);

        //- Return the word with the first letter capitalised
        inline word capitalise() const;


    // Member Operators

        // Assignment

            inline void operator=(const word&);
            inline void operator=(const string&);
            inline void operator=(const std::string&);
            inline void operator=(const char*);


    // Friend Operators

        friend word operator&(const word&, const word&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, word&);
        friend Ostream& operator<<(Ostream&, const word&);
};


void writeEntry(Ostream& os, const word& value);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "wordI.H"

#endif
