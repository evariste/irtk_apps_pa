#include "zeta.h"

Zeta::Zeta()
{
  _kZeta = 3;
  _target = NULL;
  _reference = NULL;
  _mask = NULL;
  _patchRadius = -1;
  _nbhdRadius = -1;
  _initialised = false;
  _refCount = 0;
  _patchCentreIndices = NULL;
}


Zeta::~Zeta()
{
  delete [] _reference;
  delete [] _patchCentreIndices;
  delete [] _patchOffsets;
  delete [] _nbhdOffsets;

}

void Zeta::SetTarget(irtkRealImage* tgt)
{
  _target = tgt;
}

void Zeta::SetMask(irtkGreyImage* mask)
{
  _mask = mask;
}

void Zeta::SetReferences(int count, irtkRealImage** refImg)
{
  _refCount = count;
  _reference = new irtkRealImage *[_refCount];
  for (int i = 0; i < _refCount; i++){
    _reference[i] = refImg[i];
  }
}

void Zeta::SetPatchRadius(int p)
{
  _patchRadius = p;
  cout << "Zeta::SetPatchRadius: Set patch radius to " << _patchRadius << endl;
}

void Zeta::SetNeighbourhoodRadius(int n)
{
  _nbhdRadius = n;
  cout << "Zeta::SetNeighbourhoodRadius: Set neighbourhood radius to " << _nbhdRadius << endl;
}

void Zeta::Initialise()
{
  int xdim, ydim, zdim, tdim;


  if (_target == NULL){
    cerr << "Zeta::Initialise: Target not set. Exiting." << endl;
    exit(1);
  }

  if (_refCount == 0){
    cout << "Zeta::Initialise: No reference images. Nothing to do." << endl;
    cout << "Exiting" << endl;
    exit(1);
  }

  if (_patchRadius < 0){
    cout << "Zeta::Initialise: Setting patch size to 3" << endl;
    _patchRadius = 3;
  }

  if (_nbhdRadius < 0){
    cout << "Zeta::Initialise: Setting neighbourhood size to 4" << endl;
    _nbhdRadius = 4;
  }

  xdim = _target->GetX();
  ydim = _target->GetY();
  zdim = _target->GetZ();
  tdim = _target->GetT();

  if (_mask == NULL){
    cout << "Zeta::Initialise: Create default mask" << endl;
    irtkImageAttributes attr = _target->GetImageAttributes();
    _mask = new  irtkGreyImage(attr);
    (*_mask) *= 0;
    (*_mask) += 1;
  }


  if ((xdim != _mask->GetX()) ||
      (ydim != _mask->GetY()) ||
      (zdim != _mask->GetZ()) ){
    cerr << "Zeta::Initialise: Mask dimensions don't match target. Exiting." << endl;
    exit(1);
  }


  for (int n = 0; n < _refCount; n++){
    irtkRealImage *im = _reference[n];
    if ((xdim != im->GetX()) ||
        (ydim != im->GetY()) ||
        (zdim != im->GetZ()) ||
        (tdim != im->GetT()) ){
      cerr << "Zeta::Initialise: Reference image dimensions don't match target. Exiting." << endl;
      exit(1);
    }

  }



  // The limits on where a patch centre can be.
  // Allowable indices between the low (inclusive) and the high (exclusive).
  int fovI_lo = _patchRadius + _nbhdRadius;
  int fovJ_lo = _patchRadius + _nbhdRadius;
  int fovK_lo = _patchRadius + _nbhdRadius;
  int fovI_hi = xdim - _patchRadius - _nbhdRadius;
  int fovJ_hi = ydim - _patchRadius - _nbhdRadius;
  int fovK_hi = zdim - _patchRadius - _nbhdRadius;

  if ((fovI_lo < 0) ||
      (fovI_hi > xdim) ||
      (fovJ_lo < 0) ||
      (fovJ_hi > ydim) ||
      (fovK_lo < 0) ||
      (fovK_hi > zdim))
  {
    cout << "Image to small for required neighbourhood and patch radii" << endl;
    exit(1);
  }


  // Where can we place patches?
  int n = 0;
  for (int k = fovK_lo; k < fovK_hi; k++){
    for (int j = fovJ_lo; j < fovJ_hi; j++){
      for (int i = fovI_lo; i < fovI_hi; i++){
        if ((*_mask)(i,j,k) > 0)
          n++;
      }
    }
  }

  _nPatchCentres = n;
  cout << "No of patch centres: " << _nPatchCentres << endl;

  _patchCentreIndices = new long[_nPatchCentres];

  n = 0;
  for (int k = fovK_lo; k < fovK_hi; k++){
    for (int j = fovJ_lo; j < fovJ_hi; j++){
      for (int i = fovI_lo; i < fovI_hi; i++){
        if ((*_mask)(i,j,k) > 0){
          _patchCentreIndices[n] = _target->VoxelToIndex(i, j, k);
          n++;
        }
      }
    }
  }



  // Offsets for a patch.
  cout << "Zeta::Initialise: Patch radius = " << _patchRadius << endl;

  int patchDiam = 1 + 2 * _patchRadius;
  _patchVol =  patchDiam * patchDiam * patchDiam;
  _patchOffsets = new long[_patchVol];

  // Arbitrary. Choose a centre point against which to find offsets.
  long centreIndex = _target->VoxelToIndex(fovI_lo, fovJ_lo, fovK_lo);

  n = 0;
  for (int k = -_patchRadius; k < 1 + _patchRadius; k++){
    for (int j = -_patchRadius; j < 1 + _patchRadius; j++){
      for (int i = -_patchRadius; i < 1 + _patchRadius; i++){
        long patchVoxelIndex = _target->VoxelToIndex(fovI_lo + i, fovJ_lo + j, fovK_lo + k);
        _patchOffsets[n] = patchVoxelIndex - centreIndex;
        n++;
      }
    }
  }


  // Offsets for a neighbourhood
  cout << "Zeta::Initialise: Neighbourhood radius = " << _nbhdRadius << endl;
  int nbhdDiam = 1 + 2 * _nbhdRadius;
  _nbhdVol = nbhdDiam * nbhdDiam * nbhdDiam;
  _nbhdOffsets = new long[_nbhdVol];

  n = 0;
  for (int k = -_nbhdRadius; k < 1 + _nbhdRadius; k++){
    for (int j = -_nbhdRadius; j < 1 + _nbhdRadius; j++){
      for (int i = -_nbhdRadius; i < 1 + _nbhdRadius; i++){
        long nbhdVoxelIndex = _target->VoxelToIndex(fovI_lo + i, fovJ_lo + j, fovK_lo + k);
        _nbhdOffsets[n] = nbhdVoxelIndex - centreIndex;

//        cout << "(i,j,k)" << i << " " << j << " " << k << "  ";
//        cout << n << "    " << _nbhdOffsets[n] << endl;

        n++;
      }
    }
  }

  cout << "done. " << endl;


  // TODO: Remove this test code.
  irtkRealPixel *tptr = _target->GetPointerToVoxels();
  irtkRealPixel *tgtPatchCentre, *refPatchCentre, *refPatchVox;

  (*_target) *= 0.0;

  cout << "Looping over all patches." << endl;


  for (int n = 0; n < _nPatchCentres; n++){
    tgtPatchCentre = tptr + _patchCentreIndices[n];

    for (int m = 0; m < _nbhdVol; m++){
      refPatchCentre = tgtPatchCentre + _nbhdOffsets[m];

      for (int k = 0; k < _patchVol; k++){

        *(refPatchCentre + _patchOffsets[k]) = 121;

      }
    }
  }
  _target->Write("bla.nii.gz");



  // Standardise data.

  irtkRealImage *refIm;

  gsl_matrix *X;
  X = gsl_matrix_alloc(_refCount, _nPatchCentres);


  for (n = 0; n < _refCount; n++){
    refIm = _reference[n];
    for (int k = 0; k < _nPatchCentres; k++){

    }
  }



  _initialised = true;

}











