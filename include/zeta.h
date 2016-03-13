#ifndef ZETA_H
#define ZETA_H


#include "irtkImage.h"

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

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

    void SetK(int k){_kZeta = k;}

    void Initialise();

    void Print();

    void Run();

private:

    // Number of neighbours for zeta estimation.
    int _kZeta;


    irtkRealImage *_target;

    irtkRealImage *_output;

    irtkGreyImage *_mask;

    int _refCount;
    irtkRealImage **_reference;

    int _patchRadius;

    int _nbhdRadius;

    long *_patchCentreIndices;
    int _nPatchCentres;

    long *_patchOffsets;

    long *_nbhdOffsets;

    long _tOffset;

    int _patchVol;
    int _nbhdVol;

    gsl_matrix *_Prec;

    bool _initialised;


};

#endif // ZETA_H

