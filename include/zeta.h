#ifndef ZETA_H
#define ZETA_H


#include "irtkImage.h"

#include <gsl/gsl_statistics.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>

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

    irtkRealImage *GetOutput();

private:

    void print_matrix(const gsl_matrix *m);


    // Number of neighbours for zeta estimation.
    int _kZeta;


    irtkRealImage *_target;

    irtkRealImage *_output;

    irtkGreyImage *_mask;

    int _refCount;
    irtkRealImage **_reference;

    int _patchRadius;

    int _nbhdRadius;

    unsigned long *_patchCentreIndices;
    unsigned int *_patchCentresI;
    unsigned int *_patchCentresJ;
    unsigned int *_patchCentresK;

    int _nPatchCentres;

    unsigned long *_patchOffsets;

    unsigned long *_nbhdOffsets;

    unsigned long _chanOffset;

    int _patchVol;
    int _nbhdVol;

    gsl_matrix *_Prec;

    bool _initialised;


};

#endif // ZETA_H

