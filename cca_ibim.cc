#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <fstream.h>
#include <stdio.h>
#include <string.h>
#include <irtkImage.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>
#include <nr.h>
#include <nrutil.h>

// vtk includes
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCell.h>
#include "vtkCellArray.h"
#include <vtkDataObjectWriter.h>
#include <vtkDataArray.h>
#include <vtkDataWriter.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

vtkPointSet  *data1;
vtkPoints    *points1;
vtkDataArray *vectors1;

vtkPointSet  *data2;
vtkPoints    *points2;
vtkDataArray *vectors2;

vtkPointSet  *data;
vtkPoints    *points;
vtkDataArray *vectors;

// Minimum norm of an eigenvector
// threshold for pca truncation
#define MIN_NORM 0.01 

////////////////////////////////////////////////////////////////////
// Read data into a point set.

vtkPointSet *read(char *file)
{
  int i;
  char buffer[256];
  vtkPointSet *pset;

  ifstream is(file);
  if (!is){
    cerr << "Can't open file " << file << endl;
    exit(1);
  }

  for (i = 0; i < 3; ++i){
    is.getline(buffer, 256);
  }
  is >> buffer >> buffer;

  if (strcmp(buffer, "POLYDATA") == 0){
    // Read vtkPolyData object
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(file);
    reader->Update();
    pset = reader->GetOutput();
    pset->Register(pset);
    reader->Delete();
  } else {
    if (strcmp(buffer, "UNSTRUCTURED_GRID") == 0) {
      // Read vtkUnstructuredGrid object
      vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
      reader->SetFileName(file);
      reader->Update();
      pset = reader->GetOutput();
      pset->Register(pset);
      reader->Delete();
    } else {
      if (strcmp(buffer, "STRUCTURED_GRID") == 0) {
        // Read vtkStructuredGrid object
        vtkStructuredGridReader *reader = vtkStructuredGridReader::New();
        reader->SetFileName(file);
        reader->Update();
        pset = reader->GetOutput();
        pset->Register(pset);
        reader->Delete();
      } else {
        cerr << "Unknown VTK data type" << endl;
        exit(1);
      }
    }
  }
  return pset;
}

void usage(){
  cerr << "usage: cca_ibim N";
  cerr << " struc_1_subj_1 struc_1_subj_2... struc_1_subj_N struc_2_subj_1 struc_2_subj_2 ... struc_2_subj_N <options>" << endl;
  cerr << " where N is the number of datasets." << endl;
  cerr << " Options:" << endl;
  cerr << " -modes number of modes to retain during PCA analysis of the raw data, default = 56." << endl;
  cerr << " -csv   give output in csv format:" << endl;
  cerr << "        datasets,points1,points2,pcaModes,pcaVar1,pcaVar2,pcaErr,{C|R}1,{P|N}1,C|R}2,{P|N}2,ccaCorrMax,ccaCorrAccum,ccaAccumFrac,ccaCorrAv" << endl;
  exit(1);
}

////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // this code performs a canonical correlation analysis on a number of
  // datasets.  The datapoints, i.e. the landmark data, are partitioned
  // into 2 sets, and the covariance matrix of each of these subvectors is
  // used to give a new basis.

  // If we assume that the original data is m rows (ie. m landmark
  // coordinates = 3 * no. of landmarks) by n (no of datasets), this
  // involves determining the SVD of a p by q covariance matrix where p is
  // no of landmark coords in first subvector, q is no of landmark coords
  // in second subvector (m = p + q).

  // A problem is that a p by q matrix may be very large so calculating SVD
  // may be difficult therefore we reduce the dimensions of the original
  // data first by doing a PCA on p landmarks, taking the first modesUsed
  // (modesUsed << p) modes and then rewriting the covariance data for
  // these landmarks using just these modes so that this data is now given
  // by a modesUsed by m matrix in the PCA basis.

  // Similarly we do a PCA on q landmarks and rewrite data giving another
  // modesUsed by m matrix note that vectors for m dimensional data given
  // by the composition of the individual basis vectors for p and q will be
  // orthogonal by construction, orthonormality can be achieved by
  // appropriate scaling SVD now performed on this covariance matrix
  // instead.

  int i, j, k, ok, verbose = True;
  int iNoOfLandmarks_1, iNoOfLandmarks_2;
  int iNoOfDatasets;

  // modesUsed is number of modes that remain after pca truncation, default to 56.
  int modesUsed = 56;

  irtkMatrix M, M_1, M_2;
  irtkMatrix T, T_1, T_2;
  irtkMatrix Eigenvector, Eigenvector_1, Eigenvector_2;
  irtkMatrix B, B_Eigenvector;
  irtkMatrix reducedData_1, reducedData_2;
  irtkMatrix C_11, C_12, C_21, C_22;
  irtkMatrix U, V;
  irtkMatrix Finaleigenvector, Finaleigenvector_1, Finaleigenvector_2;
  irtkMatrix D_1, D_2, D;
  irtkMatrix tempMatrix;

  irtkVector MeanShape, MeanShape_1, MeanShape_2;
  irtkVector Eigenvalues, Eigenvalues_1, Eigenvalues_2;
  irtkVector W;

  double evecNorm;
  float transformedback, check;
  double minNorm, maxNorm;
  float scaleup;
  double p[3];
  float fTotalVar   = 0;
  float explainedVariance = 0;
  char **inputNames1;
  char **inputNames2;

  if ((argc < 7)){
    usage();
  }

  // Number of datasets.
  iNoOfDatasets = atoi(argv[1]);
  argv++;
  argc--;

  inputNames1 = new char*[iNoOfDatasets];
  inputNames2 = new char*[iNoOfDatasets];

  for (i = 0; i < iNoOfDatasets; ++i){
    inputNames1[i] = argv[1];
    argv++;
    argc--;
  }

  for (i = 0; i < iNoOfDatasets; ++i){
    inputNames2[i] = argv[1];
    argv++;
    argc--;
  }

  // Read optional arguments.
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-modes") == 0)){
      argc--;
      argv++;
      modesUsed = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-csv") == 0)){
      argc--;
      argv++;
      verbose = False;
      ok = True;
    }
    if (ok == False){
      cout << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  } 

  ///////////////////////////////////////////////////////

  iNoOfLandmarks_1 = iNoOfLandmarks_2 = 0;

  // Read all landmarks_1 in matrix M_1
  for (i = 0; i < iNoOfDatasets; ++i){

    data1    = read( inputNames1[i] );

    points1  = data1->GetPoints();
    vectors1 = data1->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks_1 = points1->GetNumberOfPoints();

      M_1.Initialize(3 * iNoOfLandmarks_1, iNoOfDatasets);
      MeanShape_1.Initialize(3 * iNoOfLandmarks_1);
    }

    for (j = 0; j < iNoOfLandmarks_1; ++j){
      points1->GetPoint(j, p);

      M_1(3 * j    , i) = p[0];
      M_1(3 * j + 1, i) = p[1];
      M_1(3 * j + 2, i) = p[2];

      MeanShape_1(3 * j    ) += p[0];
      MeanShape_1(3 * j + 1) += p[1];
      MeanShape_1(3 * j + 2) += p[2];
    }
  }

 // Read all landmarks_2 in matrix M_2
  for (i = 0; i < iNoOfDatasets; ++i){

    data2    = read( inputNames2[i] );

    points2  = data2->GetPoints();
    vectors2 = data2->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks_2 = points2->GetNumberOfPoints();

      M_2.Initialize(3 * iNoOfLandmarks_2, iNoOfDatasets);
      MeanShape_2.Initialize(3 * iNoOfLandmarks_2);
    }

    for (j = 0; j < iNoOfLandmarks_2; ++j){
      points2->GetPoint(j,p);

      M_2(3 * j    , i) = p[0];
      M_2(3 * j + 1, i) = p[1];
      M_2(3 * j + 2, i) = p[2];

      MeanShape_2(3 * j    ) += p[0];
      MeanShape_2(3 * j + 1) += p[1];
      MeanShape_2(3 * j + 2) += p[2];
    }
  }

  MeanShape_1 = MeanShape_1 * (1.0 / iNoOfDatasets);
  MeanShape_2 = MeanShape_2 * (1.0 / iNoOfDatasets);

  // Subtract the mean from the data.
  for (i = 0; i < 3 * iNoOfLandmarks_1; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      M_1(i, j)  -= MeanShape_1(i);
    }
  }

  for (i = 0; i < 3 * iNoOfLandmarks_2; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      M_2(i, j)  -= MeanShape_2(i);
    }
  }

  /////////////////////////////////////////////////////////

  // Some blurb.
  if (verbose == True){
    cout << iNoOfDatasets << " datasets. " << endl;
    cout << "Landmarks :" << endl;
    cout << "Structure 1             : " << iNoOfLandmarks_1 << endl;
    cout << "Structure 2             : " << iNoOfLandmarks_2 << endl;
    cout << "Retained PCA modes      = " << modesUsed << endl;
  } else {
    cout << iNoOfDatasets << "," << iNoOfLandmarks_1 << "," << iNoOfLandmarks_2 << "," << modesUsed << "," ;
  }


  ////////////////////////////////////////////////////
  // Need to do PCA of M_1 and M_2 separately first.//
  ////////////////////////////////////////////////////

  // Form matrix 
  //   T = (1 / iNoOfDatasets) M^T * M 
  // which has dimensions iNoOfDatasets x iNoOfDatasets
  T_1.Initialize(iNoOfDatasets, iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      T_1(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_1; ++k){
        T_1(i, j) += M_1(k, i) * M_1(k, j);
      }
    }
  }
  T_1 /= iNoOfDatasets;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea_1(iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; ++i){
    for(j = 0; j < iNoOfDatasets; ++j){
      ea_1.Matrix(i, j) = T_1(i, j);
    }
  }

  ea_1.DecrSortEigenStuff();

  if (ea_1.error){
    if (verbose == True){
      cout << "ea_1.error = "<< ea_1.error << endl;
    } else {
      cout << "ea_1.err="<< ea_1.error << ",";
    }
  }

  // Convert the eigenvectors of T to true eigenvectors of the covariance
  // matrix C = M*M^T.
  //   Eigenvector_1 = M_2 * evec(T_1);
  // It is quicker to do this with the following method rather than using
  // the multiplication routines in irtkMatrix.
  Eigenvector_1.Initialize(3 * iNoOfLandmarks_1, iNoOfDatasets);

  for (i = 0; i < 3 * iNoOfLandmarks_1; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){

      Eigenvector_1(i, j) = 0.0;
      for (k = 0; k < iNoOfDatasets; ++k){
        Eigenvector_1(i, j) += M_1(i, k) * ea_1.Eigenvector(k, j);
      }
    }
  }

  // Grab the eigenvalues.
  Eigenvalues_1.Initialize(iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets ; ++i){
    Eigenvalues_1(i) = ea_1.Eigenvalue(i);
  }

  // Some reporting:
  fTotalVar   = 0;
  explainedVariance = 0;

  for (i = 0; i < iNoOfDatasets; ++i){
    fTotalVar += ea_1.Eigenvalue(i);
  }
  for (i = 0; i < modesUsed; ++i){
    explainedVariance += ea_1.Eigenvalue(i);
  }

  //////////////////////////////////////////////////////

  minNorm = FLT_MAX;
  maxNorm = -1.0 * FLT_MAX;

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets ; ++i){
    evecNorm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks_1; ++j){
      evecNorm += Eigenvector_1(j, i) * Eigenvector_1(j, i);
    }

    evecNorm = sqrt(evecNorm);

    if (evecNorm > maxNorm)
      maxNorm = evecNorm;
    if (evecNorm < minNorm)
      minNorm = evecNorm;

    // Previously only did the following if the following condition was
    // satisfied:
    //  if (100 * ea_1.Eigenvalue(i) / fTotalVar > MIN_NORM){
    // Otherwise, the eigenvector and its corresponding eigenvalue were set
    // to zero.
    for (j = 0; j < 3 * iNoOfLandmarks_1; ++j){
      Eigenvector_1(j, i) /= evecNorm;
    }
  }

  if (verbose == True){
    cout << "min and max norms 1     = " << minNorm << " " << maxNorm << endl;
    cout << "Var Struct 1            = " << 100.0 * explainedVariance / fTotalVar << "%" << endl;
  } else {
    cout << minNorm << "," << maxNorm << ",";
    cout << 100.0 * explainedVariance / fTotalVar << "%,";
  }

  ///////////////////////////////////////////////////

  // Now do the second.
  T_2.Initialize(iNoOfDatasets, iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      T_2(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_2; ++k){
        T_2(i, j) += M_2(k, i) * M_2(k, j);
      }
    }
  }
  T_2 /= iNoOfDatasets;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea_2(iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; ++i){
    for(j = 0; j < iNoOfDatasets; ++j){
      ea_2.Matrix(i, j) = T_2(i, j);
    }
  }

  ea_2.DecrSortEigenStuff();

  if (ea_2.error){
    if (verbose == True){
      cout << "ea_2.error = " << ea_2.error << endl;
    } else {
      cout << "ea_2.err=" << ea_2.error << ",";
    }
  }

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M*M^T.
  //   Eigenvector_2 = M_2 * evec(T_2);
  // It is quicker to do this with the following method rather than using
  // the multiplication routines in irtkMatrix.
  Eigenvector_2.Initialize(3 * iNoOfLandmarks_2, iNoOfDatasets);
  for (i = 0; i < 3 * iNoOfLandmarks_2; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){

      Eigenvector_2(i, j) = 0.0;
      for (k = 0; k < iNoOfDatasets; ++k){
        Eigenvector_2(i, j) += M_2(i, k) * ea_2.Eigenvector(k, j);
      }
    }
  }

  Eigenvalues_2.Initialize(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets ; ++i){
    Eigenvalues_2(i) = ea_2.Eigenvalue(i);
  }

  // Some reporting.
  fTotalVar   = 0;
  explainedVariance = 0;

  for (i = 0; i < iNoOfDatasets; ++i){
    fTotalVar += ea_2.Eigenvalue(i);
  }
  for (i = 0; i < modesUsed; ++i){
    explainedVariance += ea_2.Eigenvalue(i);
  }

  ////////////////////////////////////////////////////////

  minNorm = FLT_MAX;
  maxNorm = -1.0 * FLT_MAX;

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets ; ++i){
    evecNorm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks_2; ++j){
      evecNorm += Eigenvector_2(j, i) * Eigenvector_2(j, i);
    }

    evecNorm = sqrt(evecNorm);

    if (evecNorm > maxNorm)
      maxNorm = evecNorm;
    if (evecNorm < minNorm)
      minNorm = evecNorm;

    // Legacy code that created problems: if (100 * ea_2.Eigenvalue(i) / fTotalVar > MIN_NORM){
    if (1){
      for (j = 0; j < 3 * iNoOfLandmarks_2; ++j){
        Eigenvector_2(j, i) /= evecNorm;
      }
    } else {
      for (j = 0; j < 3 * iNoOfLandmarks_2; ++j){
        Eigenvector_2(j, i) = 0;
      }
      Eigenvalues_2(i) = 0;
    }
  }

  if (verbose == True){
    cout << "min and max norms 2     = " << minNorm << " " << maxNorm << endl;
    cout << "Var Struct 2            = " << 100.0 * explainedVariance / fTotalVar << "%" << endl;
  } else {
    cout << minNorm << "," << maxNorm << ",";
    cout << 100.0 * explainedVariance / fTotalVar << "%,";
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////


  // Now convert x-u data into pca coordinates using only first modesUsed eigenvectors.
  reducedData_1.Initialize(modesUsed, iNoOfDatasets);
  for (i = 0; i <  modesUsed; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      reducedData_1(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_1; ++k)
        reducedData_1(i, j) += Eigenvector_1(k, i) * M_1(k, j);
    }
  }

  reducedData_2.Initialize(modesUsed, iNoOfDatasets);
  for (i = 0; i <  modesUsed; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      reducedData_2(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_2; ++k)
        reducedData_2(i, j) += Eigenvector_2(k, i) * M_2(k, j);
    }
  }

  // Following are checks to see whether projections into eigenspaces worked.
  // Check data transformation is working by transforming reduced data back
  // and comparing with original data.

  if (verbose == True){
    cout << "Back projected PCA representation compared with original : " << endl;
  }

  ////// First structure:
  check = 0.0;

  for (i = 0; i < 3 * iNoOfLandmarks_1; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      transformedback = 0.0;
      for (k = 0; k < modesUsed; ++k){
        transformedback += Eigenvector_1(i, k) * reducedData_1(k, j);
      }
      check += ( transformedback - M_1(i, j) ) * ( transformedback - M_1(i, j) );
    }
  }

  if (verbose == True){
    cout << "Average squared error 1 = " << check / iNoOfLandmarks_1 / iNoOfDatasets << endl;
  } else {
    cout << check / iNoOfLandmarks_1 / iNoOfDatasets << ",";
  }

  ////// Second structure:
  check = 0.0;

  for (i = 0; i < 3 * iNoOfLandmarks_2; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      transformedback = 0.0;
      for (k = 0; k < modesUsed; ++k){
        transformedback += Eigenvector_2(i, k) * reducedData_2(k, j);
      }
      check += ( transformedback - M_2(i, j) ) * ( transformedback - M_2(i, j) );
    }
  }

  if (verbose == True){
    cout << "Average squared error 2 = " << check / iNoOfLandmarks_2 / iNoOfDatasets << endl;
  } else {
    cout << check / iNoOfLandmarks_2 / iNoOfDatasets << ",";
  }

  ///////////////////////////////////////////////////////////////////////////////////

  // Now perform CCA on the transformed, reduced dimension PCA x-u data

  // Calculate covariance matrices for reduced_data

  // Covariance blocks C_11 and C_12 can just use the eigenvalues from the
  // original data because they represent the variances along the basis
  // formed by the PCA modes.
  C_11.Initialize(modesUsed, modesUsed);
  for (i = 0; i < modesUsed; ++i){
    for (j = 0; j < modesUsed; ++j){
      C_11(i, j) = 0.0;
    }
    C_11(i, i) = ea_1.Eigenvalue(i);
  }

  C_22.Initialize(modesUsed, modesUsed);
  for (i = 0; i < modesUsed; ++i){
    for (j = 0; j < modesUsed; ++j){
      C_22(i, j) = 0.0;
    }
    C_22(i, i) = ea_2.Eigenvalue(i);
  }

  // Off-diagonal blocks need to be calculated however.
  C_12.Initialize(modesUsed, modesUsed);
  for (i = 0; i < modesUsed; ++i){
    for (j = 0; j < modesUsed; ++j){
      C_12(i, j) = 0.0;
      for (k = 0; k < iNoOfDatasets; ++k)
        C_12(i, j) += ( reducedData_1(i, k) ) * ( reducedData_2(j, k) );
    }
  }
  C_12 /= iNoOfDatasets;

  C_21.Initialize(modesUsed, modesUsed);
  for (i = 0; i < modesUsed; ++i){
    for (j = 0; j < modesUsed; ++j){
      C_21(i, j) = C_12(j, i);
    }
  }

  // perform CCA on x-u data in pca basis
  D_1.Initialize(modesUsed,modesUsed);
  D_2.Initialize(modesUsed,modesUsed);
  D.Initialize(modesUsed,modesUsed);



//    reducedData_1.Write("reducedData1.mat");
//    C_11.Write("c11.mat");
//    reducedData_2.Write("reducedData2.mat");
//    C_22.Write("c22.mat");

  if (verbose == True){
    cout << "Block determinants : " << endl;
    cout << "C_11, det = " << C_11.Det() << endl;
    cout << "C_22, det = " << C_22.Det() << endl;
  } else {
    cout << C_11.Det() << "," << C_22.Det() << ",";
  }


  // This scaling stuff is all a bit arbitrary ... try and deal with zero
  // determinants using a scale. In practice the determinants tend to be
  // very high anyway (> 10E20).
  scaleup = 1.1;

  if (C_11.Det() == 0.0) { // zero determinant, scale up
    // zero determinant, scale up.
    C_11 *= scaleup;
    C_11.Invert();
    C_11 *= scaleup;
    if (verbose == True){
      cout << "Scaled up C_11, det = " << C_11.Det() << endl;
    } else {
      cout << "sc," << C_11.Det() << ",";
    }
  } else {
    C_11.Invert();
  }

  if (C_22.Det()==0.0) {
  // zero determinant, scale up
    C_22 *= scaleup;
    C_22.Invert();
    C_22 *= scaleup;
    if (verbose == True){
      cout << "Scaled up C_22, det = " << C_22.Det() << endl;
    } else {
      cout << "sc," << C_22.Det() << ",";
    }
  } else {
    C_22.Invert();
  }

  ///////////////////////////////////////////////////////////

  // The main bit:
  
  D_1 = C_11 * C_12 * C_22 * C_21;

  // calculate eigenvectors of D_1
  float *er1, *ei1;
  float ** v1, **co, **lin;

  // Vector and matrix allocation functions in nrutil.h.
  er1 = vector(1, modesUsed);
  ei1 = vector(1, modesUsed);
  co  = matrix(1, modesUsed, 1, modesUsed);

  // Initialize real part of correlation as negative so that we know if
  // there has been a numerical recipes exit.
  for (i = 1; i <= modesUsed; ++i){
    er1[i] = -100;
  }

  // Copy D_1 matrix so can use NR routines to calculate e-values, D_1 not
  // neccesarily symmetric.
  for (i = 1; i <= modesUsed; ++i){
   for (j = 1; j <= modesUsed; ++j){
     co[i][j] = D_1(i - 1, j - 1);
   }
  }

  // Prebalance.
  balanc(co, modesUsed);
  elmhes(co, modesUsed);
  // hqr just calcs evalues
  hqr(co, modesUsed, er1, ei1); 

  // Sort eigenvalues, v1 would contain eigenvectors if we had calculated
  // them but they are just passed in for eigsorting.
  v1  = matrix(1, modesUsed, 1, modesUsed);
  eigsrt(er1, v1, modesUsed);

  //tolerance level for imaginary part to eigenvalues
  float ceigstol = 0.001; 
  int complexcorrelations = 0;
  int negativeCorrelations = 0;

  for (i = 1; i <= modesUsed; ++i){
    if ( fabs( ei1[i] ) > ceigstol)
      complexcorrelations++;

    if (er1[i] < 0)
      negativeCorrelations++;
  }

  if (verbose == True){
    if (complexcorrelations > 0)
      cout << "Structure 1: Complex correlations!" << endl;
    else 
      cout << "Structure 1: All real correlations." << endl;

    if (negativeCorrelations > 0)
      cout << "Structure 1: Negative correlations!" << endl;
    else
      cout << "Structure 1: All positive correlations." << endl;
  } else {
    if (complexcorrelations > 0)
      cout << "C,";
    else
      cout << "R,";

    if (negativeCorrelations > 0)
      cout << "N,";
    else
      cout << "P,";
  }

  // Also want to calculate average correlation
  float avcorrelation = 0.0;
  float totalCorrelation = 0.0;
  float accumulatedCorrelation = 0.0;

  for (i = 0; i < modesUsed; ++i) 
  {
    if ( er1[i+1] >= 0){
      totalCorrelation += sqrt( er1[i+1] );
      accumulatedCorrelation += totalCorrelation;
    }
  }

  avcorrelation = totalCorrelation / modesUsed;

  if (verbose == True){
   cout << "Largest Correlation     = " << sqrt(er1[1]) << endl;
   cout << "Accumulated Correlation = " << accumulatedCorrelation << endl;
   cout << "Accumulated fraction    = " << accumulatedCorrelation / (totalCorrelation * modesUsed) << endl;
   cout << "Average Correlation     = " << avcorrelation << endl;
  } else {
    cout << sqrt(er1[1]) << ",";
    cout << accumulatedCorrelation << ",";
    cout << accumulatedCorrelation / (totalCorrelation * modesUsed) << ",";
    cout << avcorrelation << endl;
  }

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cout << argv[0] << " needs to be compiled with the contrib and VTK library"
       << endl;
}
#endif


// Legacy code:



////////////////////////////////////////////////////////////////////

// // Write data out.

// void write(char *file, vtkPointSet *pset)
// {
//   if (pset->IsA("vtkStructuredGrid")){
//     vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
//     writer->SetInput((vtkStructuredGrid *)pset);
//     writer->SetFileName(file);
//     writer->SetFileTypeToBinary();
//     writer->Update();
//   } else {
//     if (pset->IsA("vtkUnstructuredGrid")){
//       vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
//       writer->SetInput((vtkUnstructuredGrid *)pset);
//       writer->SetFileName(file);
//       writer->SetFileTypeToBinary();
//       writer->Update();
//     } else {
//       if (pset->IsA("vtkPolyData")){
//         vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//         writer->SetInput((vtkPolyData *)pset);
//         writer->SetFileName(file);
//         writer->SetFileTypeToASCII();
//         writer->Write();
//       } else {
//         cout << "Unknown VTK data type" << endl;
//         exit(1);
//       }
//     }
//   }
// }


// void evex4(float **co,float *ev,float **v)
// {
//  /* calculates eigenvectors of co with eigenvalues in ev*/
//  /* eigenvectors stored in v*/
//  /* inverse iteration technique is used */

//  float **d1, **d2,*b,*bp,d,norm;
//  float epsilon;
//  int *indx,i,currentevalue,j;
//  float growthf;
//  long seed;
//  indx=ivector(1,modesUsed);
//  b=vector(1,modesUsed);
//  bp=vector(1,modesUsed);
//  epsilon=0.000001;
//  /*printf("co\n");
//  print_matrix(co,1,modesUsed,1,modesUsed);*/


//  d1=matrix(1,modesUsed,1,modesUsed);


//  for (currentevalue=1;currentevalue<=modesUsed;currentevalue++)
//  {
//   for (i=1;i<=modesUsed;++i)
//   {
//    for (j=1;j<=modesUsed;++j)
//    {
//     d1[i][j]=co[i][j];
//     if (i==j) d1[i][j]=d1[i][j]-ev[currentevalue]-0.00001*fabs(ev[currentevalue]);
//    }
//   }
//   growthf=150;
//   seed=-3;

//   /* perform LU decomposition of d1 */
//   if (ev[currentevalue]!=0){ /* non-zero eigenvalue */
//   ludcmp(d1,modesUsed,indx,&d);

//   /* now iterate to solve */
//   /* stop when mod(bk+1 - bk)<epsilon */
//   /* if growth factor norm(b1)is small, choose another random vector */

//   do
//   {
//   /*produce random normalized vector b*/
//    norm=0.0;
//    for (i=1;i<=modesUsed;++i) {
//     b[i]=gasdev(&seed);
//     norm+=b[i]*b[i];
//    }
//    norm=sqrt(norm);
//    for (i=1;i<=modesUsed;++i) b[i]=b[i]/norm;
//    lubksb(d1,modesUsed,indx,b);
//    norm=0.0;
//    for (i=1;i<=modesUsed;++i) norm+=b[i]*b[i];
//   }
//   while ((norm)<growthf);

//   cout << "out of growthf" << endl;
//   i=0;
//   float error=0;
//   do
//   {
//    /* first copy b */
//    error=0.0;
//    for (j=1;j<=modesUsed;++j) bp[j]=b[j];
//    lubksb(d1,modesUsed,indx,b);
//    norm=0.0;
//    for (j=1;j<=modesUsed;++j) norm+=b[j]*b[j];
//    norm=sqrt(norm);
//    for (j=1;j<=modesUsed;++j) b[j]=b[j]/norm;
//    for (j=1;j<=modesUsed;++j) error+=(b[j]-bp[j])*(b[j]-bp[j]);
//    ++i;
//   }
//   while ((error>epsilon)&&(i<11));
//   cout << "got eignes" << endl;
//   for (i=1;i<=modesUsed;++i) v[i][currentevalue]=b[i];
//  }
//  else
//  {
//   for (i=1;i<=modesUsed;++i) v[i][currentevalue]=0;
//  }}
// }


