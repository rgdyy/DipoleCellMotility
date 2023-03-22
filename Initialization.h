#ifndef DIPOLECELL_INITIALIZATION_H
#define DIPOLECELL_INITIALIZATION_H

#include "CommonTypes.h"

class Tissue;

struct TissueInitializer {
    static void initializeTissue(Tissue *tissue);

    static void initializeTissue2DTriangularLattice(Tissue *tissue);

    static void initializeTissue1DX(Tissue *tissue);
};

#endif // DIPOLECELL_INITIALIZATION_H