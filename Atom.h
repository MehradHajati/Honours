#ifndef ATOM_H
#define ATOM_H
#define FCC 1
#define HCP 2

#include "Vector3.h"
#include "Quaternion.h"

typedef struct{

    Vector3 position;
    int id;
    int atomType;
    int structureType;
    int selected;
    Quaternion rotationQuaternion;    

} Atom;

/// @brief Creates an atom with the given values.
/// @param position The position of the atom in space.
/// @param id An identifier for the atom.
/// @param atomType The atom's element.
/// @param structureType The structure of crystal this atom is in.
/// @param rotationQuaternion A quaternion that describes this 
/// atom's orientation.
/// @return An atom.
Atom atom_new(Vector3 position, int id, int atomType, int structureType, int selected, Quaternion rotationQuaternion);

#endif // ATOM_H