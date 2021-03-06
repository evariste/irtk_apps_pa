#ifndef ZETA_H
#define ZETA_H

// #define HAS_MPI

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

    void UseMahalanobis(bool);

    void Initialise();

    void Print();

    void Run();
#ifdef HAS_MPI
    void RunParallel();
#endif

    irtkRealImage *GetOutput();

private:

    void print_matrix(const gsl_matrix *m);

    void GetCovariance(gsl_matrix *C,   gsl_matrix *X);
    void GetPrecision(gsl_matrix *C, gsl_matrix *P);


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

    bool _initialised;

    bool _use_mahalanobis;


};

#endif // ZETA_H

