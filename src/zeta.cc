#include "zeta.h"

Zeta::Zeta()
{
  _kZeta = 3;
  _target = NULL;
  _output = NULL;
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

  irtkImageAttributes attr = _target->GetImageAttributes();

  if (_mask == NULL){
    cout << "Zeta::Initialise: Create default mask" << endl;
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



  // Set the output image!
  _output = new irtkRealImage(attr);





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

//
//  // TODO: Remove this test code.
//  irtkRealPixel *tptr = _target->GetPointerToVoxels();
//  irtkRealPixel *tgtPatchCentre, *refPatchCentre;
//
//  (*_target) *= 0.0;
//
//  cout << "Looping over all patches." << endl;
//
//
//  for (int n = 0; n < _nPatchCentres; n++){
//    tgtPatchCentre = tptr + _patchCentreIndices[n];
//
//    for (int m = 0; m < _nbhdVol; m++){
//      refPatchCentre = tgtPatchCentre + _nbhdOffsets[m];
//
//      for (int k = 0; k < _patchVol; k++){
//
//        *(refPatchCentre + _patchOffsets[k]) = 121;
//
//      }
//    }
//  }
//  _target->Write("bla.nii.gz");
//


  // Standardise data.

  irtkRealImage *refIm;

  gsl_matrix *X;
  unsigned long nDataPts = _refCount*_nPatchCentres;

  X = gsl_matrix_alloc(nDataPts, tdim);

  irtkRealPixel *refPtr;
  irtkRealPixel val;

  int tOffset = xdim * ydim * zdim;

  for (n = 0; n < _refCount; n++){

    refIm = _reference[n];
    refPtr = refIm->GetPointerToVoxels();


    for (int t = 0; t < tdim; t += tOffset){

      int i = 0;

      for (int k = 0; k < _nPatchCentres; k++){

        val = *(refPtr + _patchCentreIndices[k] + t);
        gsl_matrix_set(X, (n * _nPatchCentres) + i, t, val);

        ++i;

      }
    }
  }



  // Find the mean for each channel/modality.

  gsl_vector *meanVals = gsl_vector_alloc(tdim);

  gsl_vector_view col;
  gsl_vector_view col2;
  // Zero mean each channel/modality to get covariance

  for (int j = 0; j < X->size2; j++){

    col = gsl_matrix_column(X, j);
    double meanVal = gsl_stats_mean(col.vector.data,
                                    col.vector.stride,
                                    nDataPts);

    gsl_vector_set(meanVals, j, meanVal);

    gsl_vector_add_constant(&col.vector, -1*meanVal);

  }


  gsl_matrix *Cov = gsl_matrix_alloc(tdim, tdim);

  for (int i = 0; i < X->size2; i++){
    for (int j = 0; j < X->size2; j++){
      col = gsl_matrix_column(X, i);
      col2 = gsl_matrix_column(X, j);
      double c = gsl_stats_covariance(col.vector.data,
                                      col.vector.stride,
                                      col2.vector.data,
                                      col2.vector.stride,
                                      col.vector.size);
      cout << " c " << c << endl;
      gsl_matrix_set(Cov, i, j, c);


    }
  }




  // Set the precision matrix.
  _Prec = gsl_matrix_alloc(tdim, tdim);

  gsl_permutation * perm = gsl_permutation_alloc (tdim);
  int signum;

  gsl_linalg_LU_decomp(Cov, perm, &signum);
  gsl_linalg_LU_invert(Cov, perm, _Prec);



  _initialised = true;

}




void Zeta::Run(){

  // Loop over ROI voxels


  irtkRealPixel *tgtStartPtr, *tgtPatchCentre, *refNbhdCentre, *refStartPtr, *refPatchCentre;

  tgtStartPtr = _target->GetPointerToVoxels();


  // Store the target patches.
  gsl_matrix *T = gsl_matrix_alloc(_nPatchCentres, _patchVol);

  gsl_matrix *diff = gsl_matrix_alloc(1, _patchVol);

  gsl_matrix_view tgt_patch_vals;

  gsl_matrix * diffPrec, *diffPrecDiff;


  for (int n = 0; n < _nPatchCentres; n++){

    tgtPatchCentre = tgtStartPtr + _patchCentreIndices[n];
    for (int k = 0; k < _patchVol; k++){

      double val = *(tgtPatchCentre + _patchOffsets[k]);
      gsl_matrix_set(T, n, k, val);

    }
  }





  // For each reference
  for (int r = 0; r < _refCount; r++){

    refStartPtr = _reference[r]->GetPointerToVoxels();

    for (int n = 0; n < _nPatchCentres; n++){

      tgt_patch_vals = gsl_matrix_submatrix(T, n, 0, 1, _patchVol);

      refNbhdCentre = refStartPtr + _patchCentreIndices[n];

      // Current voxel is centre of the neighbourhood, loop
      // over reference patches in the neighbourhood.
      for (int m = 0; m < _nbhdVol; m++){

        for (int k = 0; k < _patchVol; k++){
          double val = *(refNbhdCentre + _patchOffsets[k]);
          gsl_matrix_set(diff, 0, k, val);
        }

        gsl_matrix_sub(diff, &(tgt_patch_vals.matrix));


        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
            _Prec, 0.0, diffPrec);

        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diffPrec,
            diff, 0.0, diffPrecDiff);

        cout << "diffPrecDiff size: " << diffPrecDiff->size1 << " " << diffPrecDiff->size2 << endl;
        exit(0);

      }
    }


  }




  // Find closest patch, store it and record distance
  // Find k closest patches over all references
  // Calculate gamma, mean of stored distances above
  // Calculate mean pairwise distance over k nearest patches
  // Calculate zeta


}


void Zeta::Print(){

  cout << "Patch radius: " << _patchRadius << endl;

  cout << "Neighbourhood radius: " << _nbhdRadius << endl;


  cout << "Number of voxels in ROI " << _nPatchCentres << endl;

  cout << "Number of reference images: " << _refCount << endl;

  cout << "Number of neighbours (k): " << _kZeta << endl;


  cout << "Precision matrix: " << endl;

  for (unsigned long int i = 0; i < _Prec->size1; i++){
    for (unsigned long int j = 0; j < _Prec->size2; j++){
      printf("%0.12f ", gsl_matrix_get(_Prec, i, j));
    }
    cout << endl;
  }

}




