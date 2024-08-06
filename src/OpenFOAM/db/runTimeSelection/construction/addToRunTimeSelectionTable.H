// THESE ARE MACROS FOR EASY INSERTION INTO RUN-TIME SELECTION TABLES

#ifndef addToRunTimeSelectionTable_H
#define addToRunTimeSelectionTable_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typename as the key
#define addToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add##thisType##argNames##ConstructorTo##baseType##Table_


//- Add to hash-table of functions with 'lookup' as the key
#define addNamedToRunTimeSelectionTable\
(baseType,thisType,argNames,lookup)                                            \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType::add##argNames##ConstructorToTable<thisType>                      \
        add_##lookup##_##thisType##argNames##ConstructorTo##baseType##Table_   \
        (#lookup)


//- Add to hash-table of functions with typename as the key
#define addRemovableToRunTimeSelectionTable\
(baseType,thisType,argNames)                                                   \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType::addRemovable##argNames##ConstructorToTable<thisType>             \
        addRemovable##thisType##argNames##ConstructorTo##baseType##Table_


//- Add to hash-table of functions with 'lookup' as the key
#define addRemovableNamedToRunTimeSelectionTable\
(baseType,thisType,argNames,lookup)                                            \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType::addRemovable##argNames##ConstructorToTable<thisType>             \
        addRemovable_##lookup##_##thisType##argNames##ConstructorTo            \
        ##baseType##Table_(#lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typename as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addTemplateToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames)                                              \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType::add##argNames##ConstructorToTable<thisType<Targ>>                \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Table_


//- Add to hash-table of functions with 'lookup' as the key.
//  Use when baseType doesn't need a template argument (eg, is a typedef)
#define addNamedTemplateToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookup)                                       \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType::add##argNames##ConstructorToTable<thisType<Targ>>                \
        add_##lookup##_##thisType##Targ##argNames##ConstructorTo##baseType     \
        ##Table_(#lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Add to hash-table of functions with typename as the key.
//  Use when baseType requires the Targ template argument as well
#define addTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames)                                              \
                                                                               \
    /* Add the thisType constructor function to the table */                   \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add##thisType##Targ##argNames##ConstructorTo##baseType##Targ##Table_


//- Add to hash-table of functions with 'lookup' as the key.
//  Use when baseType requires the Targ template argument as well
#define addNamedTemplatedToRunTimeSelectionTable\
(baseType,thisType,Targ,argNames,lookup)                                       \
                                                                               \
    /* Add the thisType constructor function to the table, find by lookup */   \
    baseType<Targ>::add##argNames##ConstructorToTable<thisType<Targ>>          \
        add_##lookup##_##thisType##Targ##argNames##ConstructorTo##baseType##   \
        Targ##Table_(#lookup)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
