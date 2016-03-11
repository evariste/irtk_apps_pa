#ifndef ZETA_H
#define ZETA_H


#include "irtkImage.h"

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>

class Zeta
{
public:
    Zeta();

    ~Zeta();

    void SetTarget(irtkRealImage*);

    void SetMask(irtkGreyImage *);

    void SetReferences(int, irtkRealImage**);

    void SetPatchRadius(int);

    void SetNeighbourhoodRadius(int);


    void Initialise();

private:

    // Number of neighbours for zeta estimation.
    int _kZeta;


    irtkRealImage *_target;

    irtkGreyImage *_mask;

    int _refCount;
    irtkRealImage **_reference;

    int _patchRadius;

    int _nbhdRadius;

    long *_patchCentreIndices;
    int _nPatchCentres;

    long *_patchOffsets;

    long *_nbhdOffsets;

    int _patchVol;
    int _nbhdVol;

    bool _initialised;


};

#endif // ZETA_H

