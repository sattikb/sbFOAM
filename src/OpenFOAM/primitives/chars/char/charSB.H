// CHARACTER AND POINTER TO CHARACTER STRING

#ifndef charSB_H
#define charSB_H

namespace Foam
{

// FORWARD DEFINITION OF THE FOLLOWING CLASSES.
class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
char readChar(Istream&);			// FORWARD DECLARATION OF THE FUNCTION
Istream& operator>>(Istream&, char&);
Ostream& operator<<(Ostream&, const char); 	//FOR OUTPUT OF A CHARACTER
Ostream& operator<<(Ostream&, const char*);	//FOR OUTPUT OF A STRING

inline bool isspace(char c)
{
    return
    (
        c == ' '
     || c == '\n'
     || c == '\r'
     || c == '\t'
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
