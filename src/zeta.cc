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
  _Prec = NULL;
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

  if ((fovI_lo < 0) ||
      (fovI_hi > xdim) ||
      (fovJ_lo < 0) ||
      (fovJ_hi > ydim) ||
      (fovK_lo < 0) ||
      (fovK_hi > zdim))
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


  // Standardise data.

  irtkRealImage *refIm;

  gsl_matrix *X;
  unsigned long nDataPts = _refCount*_nPatchCentres;

  X = gsl_matrix_alloc(nDataPts, nChannels);

  irtkRealPixel *refPtr;
  irtkRealPixel val;

  _chanOffset = xdim * ydim * zdim;

  for (n = 0; n < _refCount; n++){

    refIm = _reference[n];
    refPtr = refIm->GetPointerToVoxels();


    int tOff = 0;

    for (int t = 0; t < nChannels; t++){


      int i = 0;

      for (int k = 0; k < _nPatchCentres; k++){

        val = *(refPtr + _patchCentreIndices[k] + tOff);

        gsl_matrix_set(X, (n * _nPatchCentres) + i, t, val);

        ++i;

      }

      tOff += _chanOffset;


    }
  }



  // Find the mean for each channel/modality.

  gsl_vector *meanVals = gsl_vector_alloc(nChannels);

  gsl_vector_view col;
  gsl_vector_view col2;
  // Zero mean each channel/modality to get covariance

  for (unsigned int j = 0; j < X->size2; j++){

    col = gsl_matrix_column(X, j);
    double meanVal = gsl_stats_mean(col.vector.data,
                                    col.vector.stride,
                                    nDataPts);

    gsl_vector_set(meanVals, j, meanVal);

    gsl_vector_add_constant(&col.vector, -1*meanVal);

  }


  gsl_matrix *Cov = gsl_matrix_alloc(nChannels, nChannels);

  for (unsigned int i = 0; i < X->size2; i++){
    for (unsigned int j = 0; j < X->size2; j++){
      col = gsl_matrix_column(X, i);
      col2 = gsl_matrix_column(X, j);
      double c = gsl_stats_covariance(col.vector.data,
                                      col.vector.stride,
                                      col2.vector.data,
                                      col2.vector.stride,
                                      col.vector.size);
      gsl_matrix_set(Cov, i, j, c);


    }
  }




  // Set the precision matrix.
  gsl_matrix *prec = gsl_matrix_alloc(nChannels, nChannels);

  gsl_permutation * perm = gsl_permutation_alloc (nChannels);
  int signum;

  gsl_linalg_LU_decomp(Cov, perm, &signum);
  gsl_linalg_LU_invert(Cov, perm, prec);


  // Precision matrix is actually a replicated precision matrix on the block diagona.
  // There is a block for each voxel in the patch
  _Prec = gsl_matrix_alloc(_patchVol * nChannels, _patchVol * nChannels);

  gsl_matrix_set_zero(_Prec);

  for (int k = 0; k < _patchVol; k++){

    gsl_matrix_view submat;

    submat = gsl_matrix_submatrix(_Prec, k*nChannels, k*nChannels, nChannels, nChannels);

    gsl_matrix_memcpy(&(submat.matrix), prec);
  }


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


  diff = gsl_matrix_alloc(1, _patchVol*nChannels);
  gsl_vector *refPatch = gsl_vector_alloc(_patchVol*nChannels);
  //  gsl_vector *refPatchA = gsl_vector_alloc(_patchVol*tdim);
  //  gsl_vector *refPatchB = gsl_vector_alloc(_patchVol*tdim);

  gsl_matrix_view refPatchA;
  gsl_matrix_view refPatchB;

  gsl_matrix *RefData = gsl_matrix_alloc(_refCount, _patchVol*nChannels);

  double val;
  double minVal = DBL_MAX;
  double *refDistsToTgt = new double[_refCount];
  unsigned long int *sortInds = new unsigned long int[_refCount];
  double meanPairwise, meanTgtToRef;


  double normFactor = _kZeta * (_kZeta - 1) / 2.0;


  for (int i = 0; i <= 100; i++){
    printf("|");
  }
  printf("\n");
  fflush(stdout);
  int nChunk = _nPatchCentres / 100;

  // Loop over target patch centres.
  for (int n = 0; n < _nPatchCentres; n++){

    tgt_patch_vals = gsl_matrix_submatrix(T, n, 0, 1, T->size2);

    // Find closest patch to the target for each reference, store it and record distance

    // For each reference
    for (int r = 0; r < _refCount; r++){


      minVal = DBL_MAX;

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
              _Prec, 0.0, diffPrec);

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
        gsl_matrix_set_zero(diffPrec);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
            _Prec, 0.0, diffPrec);
        gsl_matrix_set_zero(dist2);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diffPrec,
            diff, 0.0, dist2);

        meanPairwise += gsl_matrix_get(dist2, 0, 0);

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


  gsl_matrix_free(diffPrec);
  gsl_matrix_free(dist2);
  gsl_matrix_free(T);
  gsl_matrix_free(diff);
  gsl_vector_free(refPatch);
  gsl_matrix_free(RefData);

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


  diff = gsl_matrix_alloc(1, _patchVol*nChannels);
  gsl_vector *refPatch = gsl_vector_alloc(_patchVol*nChannels);

  //  gsl_vector *refPatchA = gsl_vector_alloc(_patchVol*tdim);
  //  gsl_vector *refPatchB = gsl_vector_alloc(_patchVol*tdim);

  gsl_matrix_view refPatchA;
  gsl_matrix_view refPatchB;

  gsl_matrix *RefData = gsl_matrix_alloc(_refCount, _patchVol*nChannels);


  double val;
  double minVal = DBL_MAX;
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
  int myend = displs[myid] + scounts[myid];





  // Loop over target patch centres.
//  for (int n = 0; n < _nPatchCentres; n++){
  for (int n = mystart; n < myend; n++){

    tgt_patch_vals = gsl_matrix_submatrix(T, n, 0, 1, T->size2);

    // Find closest patch to the target for each reference, store it and record distance

    // For each reference
    for (int r = 0; r < _refCount; r++){


      minVal = DBL_MAX;

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
              _Prec, 0.0, diffPrec);

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
        gsl_matrix_set_zero(diffPrec);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diff,
            _Prec, 0.0, diffPrec);
        gsl_matrix_set_zero(dist2);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, diffPrec,
            diff, 0.0, dist2);

        meanPairwise += gsl_matrix_get(dist2, 0, 0);

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


  cout << "Precision matrix: " << endl;

  int nChannels = _target->GetT();

  cout << "Channels: " << nChannels << endl;

  cout << "Precision matrix " << endl;

  for (unsigned long int i = 0; i < nChannels; i++){
    for (unsigned long int j = 0; j < nChannels; j++){
      printf("%0.3f ", gsl_matrix_get(_Prec, i, j));
    }
    cout << endl;
  }

}




