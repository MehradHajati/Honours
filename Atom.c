#include "Atom.h"

Atom atom_new(Vector3 position, int id, int atomType, int structureType, int selected, Quaternion rotationQuaternion){
    Atom atom;
    atom.position = position;
    atom.id = id;
    atom.atomType = atomType;
    atom.structureType = structureType;
    atom.selected = 0;
    atom.rotationQuaternion = rotationQuaternion;
    return atom;
}