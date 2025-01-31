#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <stdio.h>
#include <string.h>

#define MAX_KEYS 100000
#define MAX_VALUE_LENGTH 100

// Structure for the inner dictionary
typedef struct {
    char key[MAX_VALUE_LENGTH];
    char value[MAX_VALUE_LENGTH];
} InnerDict;

// Structure for the outer dictionary
typedef struct {
    char key[MAX_VALUE_LENGTH];
    char value[MAX_VALUE_LENGTH];
    InnerDict innerDict[MAX_KEYS];
    int innerCount; // To track the number of inner items
} OuterDict;

#endif // DICTIONARY_H
