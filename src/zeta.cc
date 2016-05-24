#include "zeta.h"

#ifdef HAS_MPI

#include <mpi.h>

int myid;

#endif

void Zeta::print_matrix(const gsl_matrix *m)
{
    for (size_t i = 0; i < m->size1; i++) {
        for (size_t j = 0; j < m->size2; j++) {

            cout << gsl_matrix_get(m, i, j) << " ";
        }
	cout << endl;
    }
}



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
  _chanOffset = 0;
  _nPatchCentres = 0;
  _patchOffsets = NULL;
  _patchCentresI = NULL;
  _patchCentresJ = NULL;
  _patchCentresK = NULL;
  _nbhdOffsets = NULL;
  _nbhdVol = 0;
  _patchVol = 0;
  _patchCentreIndices = NULL;
  _use_mahalanobis = true;

#ifdef HAS_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

}


Zeta::~Zeta()
{

  delete [] _reference;
  delete [] _patchCentreIndices;
  delete [] _patchCentresI;
  delete [] _patchCentresJ;
  delete [] _patchCentresK;

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
#ifdef HAS_MPI
  if (myid == 0)
#endif
  {
  cout << "Zeta::SetPatchRadius: Set patch radius to " << _patchRadius << endl;
  }
}

void Zeta::SetNeighbourhoodRadius(int n)
{
  _nbhdRadius = n;
#ifdef HAS_MPI
  if (myid == 0)
#endif
  {
  cout << "Zeta::SetNeighbourhoodRadius: Set neighbourhood radius to " << _nbhdRadius << endl;
  }
}

void Zeta::UseMahalanobis(bool val){
  _use_mahalanobis = val;
#ifdef HAS_MPI
  if (myid == 0)
#endif
  {
    if (_use_mahalanobis)
      cout << "Zeta::UseMahalanobis: True" << endl;
    else
      cout << "Zeta::UseMahalanobis: False" << endl;
  }
}

irtkRealImage *Zeta::GetOutput()
{
  return _output;
}

void Zeta::GetCovariance(gsl_matrix *C,   gsl_matrix *X)
{
  // X: nObservations x nDims data matrix
  // C: Covariance matrix for X, nDims x nDims
  // BOTH MUST BE PRE-ALLOCATED.
  int dim = X->size2;
  int nDataPts = X->size1;

  gsl_vector_view colI;
  gsl_vector_view colJ;

  double val;

  for (unsigned int j = 0; j < dim; j++){

    colJ = gsl_matrix_column(X, j);
    val = gsl_stats_mean(colJ.vector.data,
                         colJ.vector.stride,
                         nDataPts);

    gsl_vector_add_constant(&colJ.vector, -1*val);
  }


  for (int i = 0; i < dim; i++){
    for (int j = i; j < dim; j++){
      colI = gsl_matrix_column(X, i);
      colJ = gsl_matrix_column(X, j);
      val = gsl_stats_covariance(colI.vector.data,
                                 colI.vector.stride,
                                 colJ.vector.data,
                                 colJ.vector.stride,
                                 colI.vector.size);
      gsl_matrix_set(C, i, j, val);
      gsl_matrix_set(C, j, i, val);
    }
  }

}

void Zeta::GetPrecision(gsl_matrix *C, gsl_matrix *P)
{
  // C: Input covariance matrix nDims x nDims
  // P: Precision matrix, inverse of C, nDims x nDims
  // BOTH MUST BE SQUARE, EQUAL SIZE AND PRE-ALLOCATED.

  int dim = C->size1;
  gsl_permutation * perm = gsl_permutation_alloc (dim);
  int signum;

  gsl_linalg_LU_decomp(C, perm, &signum);

  double det = gsl_linalg_LU_det(C, signum);

  if (det > 0.0001)
    gsl_linalg_LU_invert(C, perm, P);
  else
    gsl_matrix_set_identity(P);

  gsl_permutation_free(perm);
}

void Zeta::Initialise()
{
  int xdim, ydim, zdim, nChannels;

// #ifdef HAS_MPI
//   int myid;
//   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
// #endif

  if (_target == NULL){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cerr << "Zeta::Initialise: Target not set. Exiting." << endl;
    }
    exit(1);
  }

  if (_refCount == 0){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cerr << "Zeta::Initialise: No reference images. Nothing to do." << endl;
    cerr << "Exiting" << endl;
    }
    exit(1);
  }

  if (_patchRadius < 0){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cout << "Zeta::Initialise: Setting patch size to 3" << endl;
    }
    _patchRadius = 3;
  }

  if (_nbhdRadius < 0){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cout << "Zeta::Initialise: Setting neighbourhood size to 4" << endl;
    }
    _nbhdRadius = 4;
  }

  if (_kZeta < 2){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cerr << "Zeta::Initialise: Number of nearest neighbours is less than 2. Exiting" << endl;
    }
    exit(1);
  }

  xdim = _target->GetX();
  ydim = _target->GetY();
  zdim = _target->GetZ();
  nChannels = _target->GetT();

  irtkImageAttributes attr = _target->GetImageAttributes();

  if (_mask == NULL){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cout << "Zeta::Initialise: Create default mask" << endl;
    }
    _mask = new  irtkGreyImage(attr);
    (*_mask) *= 0;
    (*_mask) += 1;
  }


  if ((xdim != _mask->GetX()) ||
      (ydim != _mask->GetY()) ||
      (zdim != _mask->GetZ()) ){
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cerr << "Zeta::Initialise: Mask dimensions don't match target. Exiting." << endl;
    }
    exit(1);
  }


  // Set the output image! it only needs to be a 3D volume.
  attr._t = 1;
  _output = new irtkRealImage;
  _output->Initialize(attr);

  for (int n = 0; n < _refCount; n++){
    irtkRealImage *im = _reference[n];
    if ((xdim != im->GetX()) ||
        (ydim != im->GetY()) ||
        (zdim != im->GetZ()) ||
        (nChannels != im->GetT()) ){
#ifdef HAS_MPI
    if (myid == 0)
#endif
      {
      cerr << "Zeta::Initialise: Reference image dimensions don't match target. Exiting." << endl;
      }
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

  if ((fovI_lo < 0)        ||
      (fovI_hi > xdim)     ||
      (fovJ_lo < 0)        ||
      (fovJ_hi > ydim)     ||
      (fovK_lo < 0)        ||
      (fovK_hi > zdim)     ||
      (fovI_lo >= fovI_hi) ||
      (fovJ_lo >= fovJ_hi) ||
      (fovK_lo >= fovK_hi) )
  {
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cout << "Image to small for required neighbourhood and patch radii" << endl;
    }
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
#ifdef HAS_MPI
    if (myid == 0)
#endif
  {
  cout << "No of patch centres: " << _nPatchCentres << endl;
  }

  _patchCentreIndices = new unsigned long[_nPatchCentres];
  _patchCentresI = new unsigned int[_nPatchCentres];
  _patchCentresJ = new unsigned int[_nPatchCentres];
  _patchCentresK = new unsigned int[_nPatchCentres];

  n = 0;
  for (int k = fovK_lo; k < fovK_hi; k++){
    for (int j = fovJ_lo; j < fovJ_hi; j++){
      for (int i = fovI_lo; i < fovI_hi; i++){
        if ((*_mask)(i,j,k) > 0){
          _patchCentreIndices[n] = _target->VoxelToIndex(i, j, k);
          _patchCentresI[n] = i;
          _patchCentresJ[n] = j;
          _patchCentresK[n] = k;
          n++;
        }
      }
    }
  }



  // Offsets for a patch.
#ifdef HAS_MPI
    if (myid == 0)
#endif
  {
  cout << "Zeta::Initialise: Patch radius = " << _patchRadius << endl;
  }

  int patchDiam = 1 + 2 * _patchRadius;
  _patchVol =  patchDiam * patchDiam * patchDiam;
  _patchOffsets = new unsigned long[_patchVol];

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
#ifdef HAS_MPI
    if (myid == 0)
#endif
  {
  cout << "Zeta::Initialise: Neighbourhood radius = " << _nbhdRadius << endl;
  }
  int nbhdDiam = 1 + 2 * _nbhdRadius;
  _nbhdVol = nbhdDiam * nbhdDiam * nbhdDiam;
  _nbhdOffsets = new unsigned long[_nbhdVol];

  n = 0;
  for (int k = -_nbhdRadius; k < 1 + _nbhdRadius; k++){
    for (int j = -_nbhdRadius; j < 1 + _nbhdRadius; j++){
      for (int i = -_nbhdRadius; i < 1 + _nbhdRadius; i++){
        long nbhdVoxelIndex = _target->VoxelToIndex(fovI_lo + i, fovJ_lo + j, fovK_lo + k);
        _nbhdOffsets[n] = nbhdVoxelIndex - centreIndex;

        n++;
      }
    }
  }

#ifdef HAS_MPI
    if (myid == 0)
#endif
  {
  cout << "done. " << endl;
  }

  // Offset for moving from a voxel in one 3D volume to
  // the corresponding voxel in the next volume (assuming a 4D image)
  _chanOffset = xdim * ydim * zdim;


  _initialised = true;

#ifdef HAS_MPI
    if (myid == 0)
#endif
  {
  cout << "Done initialising" << endl;
  }

}




void Zeta::Run(){

  if (! _initialised){
    this->Initialise();
  }

  // Loop over ROI voxels


  irtkRealPixel *tgtStartPtr, *tgtPatchCentre, *refNbhdCentre, *refPatchCentre, *pTemp;

  irtkRealPixel **refStartPtr;

  refStartPtr = new irtkRealPixel*[_refCount];

  for (int r = 0; r < _refCount; r++){
    refStartPtr[r] = _reference[r]->GetPointerToVoxels();
  }



  tgtStartPtr = _target->GetPointerToVoxels();


  int nChannels = _target->GetT();


  gsl_matrix_view tgt_patch_vals;
  gsl_matrix *diff, *diffPrec, *dist2;


  diffPrec = gsl_matrix_alloc(1, _patchVol * nChannels);
  dist2 = gsl_matrix_alloc(1, 1);


  // Store the target patches.

  // Each row is of length nChannels x patch volume. A row
  // contains the values for a patch in the first channel
  // followed by the patch values for the second channel and so on.

  gsl_matrix *T = gsl_matrix_alloc(_nPatchCentres, _patchVol*nChannels);

  for (int n = 0; n < _nPatchCentres; n++){

    tgtPatchCentre = tgtStartPtr + _patchCentreIndices[n];

    for (int k = 0; k < _patchVol; k++){

      int offset = 0;

      pTemp = tgtPatchCentre + _patchOffsets[k];

      for (int chan = 0; chan < nChannels; chan++, offset += _chanOffset){

        double val = *(pTemp + offset ) ;

        gsl_matrix_set(T, n, k + chan*_patchVol, val);
      }
    }
  }

  // TODO: Transport to MPI section
  // Find per-voxel precision matrix from Reference images.
  // Only need to store upper triangle of precision matrix.

  // There will be _refCount*_patchVol samples used to estimate a
  // square covariance and precision matrix of size nChannels. Need to
  // Check we have enough samples otherwise the covariance will be singular.

  if (_refCount * _patchVol < nChannels){
    cerr << "Insufficient samples to estimate covariance. " << endl;
    cerr << "Consider increasing patch radius or number of reference images." << endl;
    exit(1);
  }

  unsigned long nPrecMatEntries = nChannels * (nChannels + 1) / 2;
  gsl_matrix* precMatData = gsl_matrix_alloc(_nPatchCentres, nPrecMatEntries);

  // TODO: Make conditional if using mahalanobis

  // Voxel data from all reference images at a patch location.
  gsl_matrix* voxData = gsl_matrix_alloc(_refCount*_patchVol, nChannels);

  gsl_vector *refVoxel = gsl_vector_alloc(nChannels);

  // Covariance:
  gsl_matrix *C = gsl_matrix_alloc(nChannels, nChannels);
  // Precision:
  gsl_matrix *P = gsl_matrix_alloc(nChannels, nChannels);

  // Loop over patch centres.
  for (int n = 0; n < _nPatchCentres; n++){

    // Get all reference patches for current patch centre.
    for (int r = 0; r < _refCount; r++){

      refPatchCentre = refStartPtr[r] + _patchCentreIndices[n];

      // Data for current reference.
      for (int k = 0; k < _patchVol; k++){

        // Particular voxel within the patch.
        pTemp = refPatchCentre + _patchOffsets[k];

        // Loop over channels
        int offset = 0;
        for (int chan = 0; chan < nChannels; chan++, offset += _chanOffset){

          gsl_vector_set(refVoxel,
                         chan,
                         *(pTemp + offset));
        }

        // Store voxel info.
        gsl_matrix_set_row(voxData,
                           r * _patchVol + k,
                           refVoxel);

      }
    } // References

    // Get precision for current voxel.
    GetCovariance(C, voxData);
    GetPrecision(C, P);

    int k = 0;
    for (int i=0; i < nChannels; i++){
      for (int j = i; j < nChannels; j++){
        double val = gsl_matrix_get(P, i, j);
        gsl_matrix_set(precMatData, n, k, val);
        k++;
      }
    }

  } // Patch centres



  // Variables for distance calculation

  gsl_vector *refPatch = gsl_vector_alloc(_patchVol*nChannels);
  gsl_matrix *RefData = gsl_matrix_alloc(_refCount, _patchVol*nChannels);

  diff = gsl_matrix_alloc(1, _patchVol*nChannels);

  gsl_matrix_view refPatchA;
  gsl_matrix_view refPatchB;


  double val, minVal;
  double *refDistsToTgt = new double[_refCount];
  unsigned long int *sortInds = new unsigned long int[_refCount];
  double meanPairwise, meanTgtToRef;

  double normFactor = _kZeta * (_kZeta - 1) / 2.0;


  // Progress markers
  for (int i = 0; i <= 100; i++){
    printf("|");
  }
  printf("\n");
  fflush(stdout);
  int nChunk = _nPatchCentres / 100;

  gsl_matrix *prec = gsl_matrix_alloc(_patchVol * nChannels, _patchVol * nChannels);


  // Loop over target patch centres.
  for (int n = 0; n < _nPatchCentres; n++){

    tgt_patch_vals = gsl_matrix_submatrix(T, n, 0, 1, T->size2);

    // Get precision matrix for current patch.
    // TODO: Make conditional if using mahalanobis

    int k = 0;
    for (int i=0; i < nChannels; i++){
      for (int j = i; j < nChannels; j++){
        double val = gsl_matrix_get(precMatData, n, k);
        gsl_matrix_set(P, i, j, val);
        gsl_matrix_set(P, j, i, val);
        k++;
      }
    }

    // Replicate the precision matrix along the diagonal of block matrix for the entire patch.
    gsl_matrix_set_zero(prec);
    for (int k = 0; k < _patchVol; k++){
      gsl_matrix_view submat;
      submat = gsl_matrix_submatrix(prec, k*nChannels, k*nChannels, nChannels, nChannels);
      gsl_matrix_memcpy(&(submat.matrix), P);
    }


    // Find closest patch to the target for each reference, store it and record distance

    // For each reference
    for (int r = 0; r < _refCount; r++){


      minVal = 100000000.0;

      refNbhdCentre = refStartPtr[r] + _patchCentreIndices[n];

      // Current voxel is centre of the neighbourhood, loop
      // over reference patches in the neighbourhood.
      for (int m = 0; m < _nbhdVol; m++){

        refPatchCentre = refNbhdCentre + _nbhdOffsets[m];

        // Data for current reference patch.
        for (int k = 0; k < _patchVol; k++){

          int offset = 0;
          // Particular spatial location within the patch.
          pTemp = refPatchCentre + _patchOffsets[k];

          for (int chan = 0; chan < nChannels; chan++, offset += _chanOffset){

            double val = *(pTemp + offset);

            gsl_vector_set(refPatch, k+chan*_patchVol, val);
          }
        }

        gsl_matrix_set_row(diff, 0, refPatch);

        gsl_matrix_sub(diff, &(tgt_patch_vals.matrix));

        if (_use_mahalanobis)
        {
          gsl_matrix_set_zero(diffPrec);

          // diff is a row vector, 1 x (patch vol * channels)
          // diffPrec = diff * _Prec
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
              prec, 0.0, diffPrec);

          gsl_matrix_set_zero(dist2);

          // dist2 = diffPrec * diff^T
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diffPrec,
              diff, 0.0, dist2);
        }
        else
        {
          // dist2 = diff * diff^T
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diff,
              diff, 0.0, dist2);
        }

        val = sqrt( gsl_matrix_get(dist2, 0, 0) );

        if (minVal > val){
          // Pick the current patch as the best one from this reference.
          minVal = val;
          gsl_matrix_set_row(RefData, r, refPatch);
        }

      } // Loop over patches in neighbourhood for current reference image.


      // Set the minimum value for the current reference.
      refDistsToTgt[r] = minVal;

    } // Loop over references


    // Gamma: Mean 'distance' to k-nearest refs.

    // Sort and take mean of first k
    gsl_sort_smallest_index(sortInds, _kZeta, refDistsToTgt, 1, _refCount);

    meanTgtToRef = 0.0;
    for (int r = 0; r < _kZeta; r++){
      meanTgtToRef += refDistsToTgt[sortInds[r]];
    }
    meanTgtToRef /= _kZeta;


    // Get pairwise distances in the k-nearest clique of references and find their mean.
    meanPairwise = 0.0;
    for (int rA = 0; rA < _kZeta-1; rA++){

      refPatchA = gsl_matrix_submatrix(RefData,
                      sortInds[rA], 0, 1, RefData->size2);

      for (int rB = rA+1; rB < _kZeta; rB++){

        refPatchB = gsl_matrix_submatrix(RefData,
                        sortInds[rB], 0, 1, RefData->size2);

        gsl_matrix_memcpy(diff, &(refPatchB.matrix));
        gsl_matrix_sub(diff, &(refPatchA.matrix));

        if (_use_mahalanobis)
        {
          // dist2 = diff * diffPrec * diff^T
          gsl_matrix_set_zero(diffPrec);
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
              prec, 0.0, diffPrec);
          gsl_matrix_set_zero(dist2);
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diffPrec,
              diff, 0.0, dist2);
        }
        else
        {
          // dist2 = diff * diff^T
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diff,
              diff, 0.0, dist2);
        }

        meanPairwise += sqrt( gsl_matrix_get(dist2, 0, 0) );

      }
    }

    meanPairwise /= normFactor;

    // Zeta is difference between Gamma and mean intra-clique distance.
    _output->PutAsDouble(_patchCentresI[n], _patchCentresJ[n], _patchCentresK[n], meanTgtToRef - meanPairwise);


    if (n % nChunk == 0){
      printf(".");
      fflush(stdout);
    }

  } // Loop over patch centres, index: n

  cout << endl << endl;

  gsl_matrix_free(diffPrec);
  gsl_matrix_free(dist2);
  gsl_matrix_free(T);
  gsl_matrix_free(diff);
  gsl_vector_free(refPatch);
  gsl_matrix_free(RefData);
  gsl_matrix_free(precMatData);

  gsl_matrix_free(voxData);
  gsl_vector_free(refVoxel);
  gsl_matrix_free(C);
  gsl_matrix_free(P);

  delete [] refStartPtr;
  delete [] refDistsToTgt;
  delete [] sortInds;
}




#ifdef HAS_MPI

void Zeta::RunParallel()
{

  if (! _initialised){
    this->Initialise();
  }


  irtkRealPixel *tgtStartPtr, *tgtPatchCentre, *refNbhdCentre, *refPatchCentre, *pTemp;

  irtkRealPixel **refStartPtr;

  refStartPtr = new irtkRealPixel*[_refCount];

  for (int r = 0; r < _refCount; r++){
    refStartPtr[r] = _reference[r]->GetPointerToVoxels();
  }

  /* get myid and # of processors */
  int myid,numproc;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (myid == 0)
    cout << "Zeta: Running parallel version." << endl;


  tgtStartPtr = _target->GetPointerToVoxels();


  int nChannels = _target->GetT();


  gsl_matrix_view tgt_patch_vals;
  gsl_matrix *diff, *diffPrec, *dist2;


  diffPrec = gsl_matrix_alloc(1, _patchVol * nChannels);
  dist2 = gsl_matrix_alloc(1, 1);





  diff = gsl_matrix_alloc(1, _patchVol*nChannels);
  gsl_vector *refPatch = gsl_vector_alloc(_patchVol*nChannels);

  //  gsl_vector *refPatchA = gsl_vector_alloc(_patchVol*tdim);
  //  gsl_vector *refPatchB = gsl_vector_alloc(_patchVol*tdim);

  gsl_matrix_view refPatchA;
  gsl_matrix_view refPatchB;

  gsl_matrix *RefData = gsl_matrix_alloc(_refCount, _patchVol*nChannels);


  double val, minVal;
  double *refDistsToTgt = new double[_refCount];
  unsigned long int *sortInds = new unsigned long int[_refCount];
  double meanPairwise, meanTgtToRef;


  double normFactor = _kZeta * (_kZeta - 1) / 2.0;



  // Array to populate data from all processes.
  double *outputValues;



  if (myid == 0)
    outputValues = new double[_nPatchCentres];


  /* divide loop */
  int *displs = new int[numproc];
  int *scounts = new int[numproc];

  int currDisp = 0;
  int chunkSizeA = _nPatchCentres / numproc;
  int chunkSizeB = 1 + chunkSizeA;

  for (int i = 0; i < numproc; i++){
    displs[i] = currDisp;

    if ( i < _nPatchCentres % numproc ){
      currDisp += chunkSizeB;
      scounts[i] = chunkSizeB;
    } else {
      currDisp += chunkSizeA;
      scounts[i] = chunkSizeA;
    }
  }

  if (myid == 0){
    for (int i = 0; i < numproc; i++){
      printf("  i: %02d  displ: %d   scount: %d\n",i, displs[i], scounts[i]);
    }
  }


  double *rbuf = new double[scounts[myid]];

  MPI_Scatterv(outputValues, scounts, displs, MPI_DOUBLE, rbuf, scounts[myid], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Compute for the following index limits */
  int mystart = displs[myid];
  int mycount = scounts[myid];
  int myend = displs[myid] + mycount;



  // Store the target patches.

  // Each row is of length nChannels x patch volume. A row
  // contains the values for a patch in the first channel
  // followed by the patch values for the second channel and so on.

  gsl_matrix *T = gsl_matrix_alloc(mycount, _patchVol*nChannels);

  // TODO: Make the range of index values restricted to those for the current process.
  for (int n = mystart; n < myend; n++){

    tgtPatchCentre = tgtStartPtr + _patchCentreIndices[n];

    for (int k = 0; k < _patchVol; k++){

      int offset = 0;

      pTemp = tgtPatchCentre + _patchOffsets[k];

      for (int chan = 0; chan < nChannels; chan++, offset += _chanOffset){

        double val = *(pTemp + offset ) ;

        gsl_matrix_set(T, n - mystart, k + chan*_patchVol, val);
      }
    }
  }




  // Find per-voxel precision matrix from Reference images.
  // Only need to store upper triangle of precision matrix.

  // There will be _refCount*_patchVol samples used to estimate a
  // square covariance and precision matrix of size nChannels. Need to
  // Check we have enough samples otherwise the covariance will be singular.

  if ((myid == 0) && (_refCount * _patchVol < nChannels)){
    cerr << "Insufficient samples to estimate covariance. " << endl;
    cerr << "Consider increasing patch radius or number of reference images." << endl;
    exit(1);
  }

  unsigned long nPrecMatEntries = nChannels * (nChannels + 1) / 2;
  gsl_matrix* precMatData = gsl_matrix_alloc(mycount, nPrecMatEntries);

  // TODO: Make conditional if using mahalanobis

  // Voxel data from all reference images at a patch location.
  gsl_matrix* voxData = gsl_matrix_alloc(_refCount*_patchVol, nChannels);

  gsl_vector *refVoxel = gsl_vector_alloc(nChannels);

  // Covariance:
  gsl_matrix *C = gsl_matrix_alloc(nChannels, nChannels);
  // Precision:
  gsl_matrix *P = gsl_matrix_alloc(nChannels, nChannels);

  // Loop over patch centres.
  for (int n = mystart; n < myend; n++){

    // Get all reference patches for current patch centre.
    for (int r = 0; r < _refCount; r++){

      refPatchCentre = refStartPtr[r] + _patchCentreIndices[n];

      // Data for current reference.
      for (int k = 0; k < _patchVol; k++){

        // Particular voxel within the patch.
        pTemp = refPatchCentre + _patchOffsets[k];

        // Loop over channels
        int offset = 0;
        for (int chan = 0; chan < nChannels; chan++, offset += _chanOffset){

          gsl_vector_set(refVoxel,
                         chan,
                         *(pTemp + offset));
        }

        // Store voxel info.
        gsl_matrix_set_row(voxData,
                           r * _patchVol + k,
                           refVoxel);

      }
    } // References

    // Get precision for current voxel.
    GetCovariance(C, voxData);
    GetPrecision(C, P);

    int k = 0;
    for (int i=0; i < nChannels; i++){
      for (int j = i; j < nChannels; j++){
        double val = gsl_matrix_get(P, i, j);
        gsl_matrix_set(precMatData, n - mystart, k, val);
        k++;
      }
    }

  } // Patch centres



  gsl_matrix *prec = gsl_matrix_alloc(_patchVol * nChannels, _patchVol * nChannels);


  // Now do the actual work. Loop over target patch centres.
  for (int n = mystart; n < myend; n++){

    tgt_patch_vals = gsl_matrix_submatrix(T, n - mystart, 0, 1, T->size2);

    // Find closest patch to the target for each reference, store it and record distance

    // For each reference
    for (int r = 0; r < _refCount; r++){


      minVal = 100000000.0;

      refNbhdCentre = refStartPtr[r] + _patchCentreIndices[n];

      // Current voxel is centre of the neighbourhood, loop
      // over reference patches in the neighbourhood.
      for (int m = 0; m < _nbhdVol; m++){

        refPatchCentre = refNbhdCentre + _nbhdOffsets[m];

        // Data for current reference patch.
        for (int k = 0; k < _patchVol; k++){

          int offset = 0;
          // Particular spatial location within the patch.
          pTemp = refPatchCentre + _patchOffsets[k];

          for (int chan = 0; chan < nChannels; chan++, offset += _chanOffset){

            double val = *(pTemp + offset);

            gsl_vector_set(refPatch, k+chan*_patchVol, val);
          }
        }

        gsl_matrix_set_row(diff, 0, refPatch);

        gsl_matrix_sub(diff, &(tgt_patch_vals.matrix));

        if (_use_mahalanobis)
        {

          // Get precision matrix for current patch.
          int k = 0;
          for (int i=0; i < nChannels; i++){
            for (int j = i; j < nChannels; j++){
              double val = gsl_matrix_get(precMatData, n - mystart, k);
              gsl_matrix_set(P, i, j, val);
              gsl_matrix_set(P, j, i, val);
              k++;
            }
          }

          // Replicate the precision matrix along the diagonal of block matrix for the entire patch.
          gsl_matrix_set_zero(prec);
          for (int k = 0; k < _patchVol; k++){
            gsl_matrix_view submat;
            submat = gsl_matrix_submatrix(prec, k*nChannels, k*nChannels, nChannels, nChannels);
            gsl_matrix_memcpy(&(submat.matrix), P);
          }


          gsl_matrix_set_zero(diffPrec);

          // diff is a row vector, 1 x (patch vol * channels)
          // diffPrec = diff * _Prec
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
              prec, 0.0, diffPrec);

          gsl_matrix_set_zero(dist2);

          // dist2 = diffPrec * diff^T
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diffPrec,
              diff, 0.0, dist2);
        }
        else
        {
          // dist2 = diff * diff^T
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diff,
              diff, 0.0, dist2);
        }

        val = sqrt( gsl_matrix_get(dist2, 0, 0) );

        if (minVal > val){
          // Pick the current patch as the best one from this reference.
          minVal = val;
          gsl_matrix_set_row(RefData, r, refPatch);
        }

      } // Loop over patches in neighbourhood for current reference image.


      // Set the minimum value for the current reference.
      refDistsToTgt[r] = minVal;

    } // Loop over references


    // Gamma: Mean 'distance' to k-nearest refs.

    // Sort and take mean of first k
    gsl_sort_smallest_index(sortInds, _kZeta, refDistsToTgt, 1, _refCount);

    meanTgtToRef = 0.0;
    for (int r = 0; r < _kZeta; r++){
      meanTgtToRef += refDistsToTgt[sortInds[r]];
    }
    meanTgtToRef /= _kZeta;


    // Get pairwise distances in the k-nearest clique of references and find their mean.
    meanPairwise = 0.0;
    for (int rA = 0; rA < _kZeta-1; rA++){

      refPatchA = gsl_matrix_submatrix(RefData,
                      sortInds[rA], 0, 1, RefData->size2);

      for (int rB = rA+1; rB < _kZeta; rB++){

        refPatchB = gsl_matrix_submatrix(RefData,
                        sortInds[rB], 0, 1, RefData->size2);

        gsl_matrix_memcpy(diff, &(refPatchB.matrix));
        gsl_matrix_sub(diff, &(refPatchA.matrix));

        if (_use_mahalanobis)
        {
          // dist2 = diff * diffPrec * diff^T
          gsl_matrix_set_zero(diffPrec);
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
              prec, 0.0, diffPrec);
          gsl_matrix_set_zero(dist2);
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diffPrec,
              diff, 0.0, dist2);
        }
        else
        {
          // dist2 = diff * diff^T
          gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diff,
              diff, 0.0, dist2);
        }

        meanPairwise += sqrt( gsl_matrix_get(dist2, 0, 0) );

      }
    }

    meanPairwise /= normFactor;

    // Zeta is difference between Gamma and mean intra-clique distance.
    rbuf[n - mystart] = meanTgtToRef - meanPairwise;


  } // Loop over patch centres, index: n


  MPI_Gatherv(rbuf, scounts[myid],
              MPI_DOUBLE, outputValues, scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  // Loop over target patch centres.
  if (myid == 0){
    for (int n = 0; n < _nPatchCentres; n++){
      _output->PutAsDouble(_patchCentresI[n], _patchCentresJ[n], _patchCentresK[n], outputValues[n]);
    }
  }


  gsl_matrix_free(diffPrec);
  gsl_matrix_free(dist2);
  gsl_matrix_free(T);
  gsl_matrix_free(diff);
  gsl_vector_free(refPatch);
  gsl_matrix_free(RefData);

  delete [] refStartPtr;
  delete [] refDistsToTgt;
  delete [] sortInds;

  if (myid == 0)
    delete [] outputValues;

}

#endif




void Zeta::Print(){

  cout << "Patch radius: " << _patchRadius << endl;

  cout << "Neighbourhood radius: " << _nbhdRadius << endl;


  cout << "Number of voxels in ROI " << _nPatchCentres << endl;

  cout << "Number of reference images: " << _refCount << endl;

  cout << "Number of neighbours (k): " << _kZeta << endl;


  int nChannels = _target->GetT();

  cout << "Channels: " << nChannels << endl;


}





