#include "Atom.h"


/**
 * @brief Creates an atom with the given values.
 * 
 * @param position The position of the atom in space.
 * @param id An identifier for the atom.
 * @param atomType The atom's element.
 * @param structureType The structure of crystal this atom is in.
 * @param selected Not used currently in this version, this number is set to zero in the struct always
 * @param rotationQuaternion A quaternion that describes this atom's orientation.
 * @return Atom struct holding the information
 */
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