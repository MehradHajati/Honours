#ifndef SAMPLESTRETCHER_H
#define SAMPLESTRETCHER_H

#include "Tilt.h"

/// @brief Stretches the given smaple's tilted data to fill the area occupied by the untilted data.
/// @param smpl The sample to stretch.
void sampleStretcher_stretch(Sample *smpl, double phi);

#endif // SAMPLESTRETCHER_H
