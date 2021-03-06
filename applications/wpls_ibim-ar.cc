#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <fstream.h>
#include <stdio.h>
#include <string.h>
#include <irtkImage.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>
#include <nr.h>
#include <nrutil.h>

// The number of modes that are kept from the pca this is based on the number of datasets;
// There are 93 subjects, leaving one out gives 92.  Could choose fewer modes if we want.


// #define noOfModes 177

// Minimum norm of an eigenvector
#define MIN_NORM 0.01

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
vtkPointSet *data1;
vtkPoints *points1;
vtkDataArray *vectors1;
vtkPointSet *data2;
vtkPoints *points2;
vtkDataArray *vectors2;


vtkPointSet *
read (char *file)
{
  int i;
  char buffer[256];
  vtkPointSet *pset;

  ifstream is (file);
  if (!is) {
    cerr << "Can't open file " << file << endl;
    exit (1);
  }

  for (i = 0; i < 3; i++) {
    is.getline (buffer, 256);
  }
  is >> buffer >> buffer;

  if (strcmp (buffer, "POLYDATA") == 0) {

    // Read vtkPolyData object
    vtkPolyDataReader *reader = vtkPolyDataReader::New ();
    reader->SetFileName (file);
    reader->Update ();
    pset = reader->GetOutput ();
    pset->Register (pset);
    reader->Delete ();
  }
  else {
    if (strcmp (buffer, "UNSTRUCTURED_GRID") == 0) {
      // Read vtkUnstructuredGrid object
      vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New ();
      reader->SetFileName (file);
      reader->Update ();
      pset = reader->GetOutput ();
      pset->Register (pset);
      reader->Delete ();
    }
    else {
      if (strcmp (buffer, "STRUCTURED_GRID") == 0) {
	// Read vtkStructuredGrid object
	vtkStructuredGridReader *reader = vtkStructuredGridReader::New ();
	reader->SetFileName (file);
	reader->Update ();
	pset = reader->GetOutput ();
	pset->Register (pset);
	reader->Delete ();
      }
      else {
	cerr << "Unknown VTK data type" << endl;
	exit (1);
      }
    }
  }
  return pset;
}

void
write (char *file, vtkPointSet * pset)
{
  if (pset->IsA ("vtkStructuredGrid")) {
    vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New ();
    writer->SetInput ((vtkStructuredGrid *) pset);
    writer->SetFileName (file);
    writer->SetFileTypeToBinary ();
    writer->Update ();
  }
  else {
    if (pset->IsA ("vtkUnstructuredGrid")) {
      vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New ();
      writer->SetInput ((vtkUnstructuredGrid *) pset);
      writer->SetFileName (file);
      writer->SetFileTypeToBinary ();
      writer->Update ();
    }
    else {
      if (pset->IsA ("vtkPolyData")) {
	vtkPolyDataWriter *writer = vtkPolyDataWriter::New ();
	writer->SetInput ((vtkPolyData *) pset);
	writer->SetFileName (file);
	writer->SetFileTypeToASCII ();
	writer->Write ();
      }
      else {
	cerr << "Unknown VTK data type" << endl;
	exit (1);
      }
    }
  }
}

float
sum_abs_matrix (irtkMatrix M)
{
  int i, j;
  float result = 0.0;
  for (i = 0; i < M.Rows (); i++) {
    for (j = 0; j < M.Cols (); j++) {
      result += fabs (M (i, j));
    }
  }
  return result;
}

void
print_matrix (irtkMatrix M)
{
  int i, j;
  for (i = 0; i < M.Rows (); i++) {
    for (j = 0; j < M.Cols (); j++) {
      cerr << M (i, j);
      if (j == M.Cols () - 1)
	cerr << endl;
      else
	cerr << " ";
    }
  }
}


void
write_matrix (irtkMatrix M)
{
  ofstream file2 ("matrixerror", ios::out | ios::app);
  //cerr << "no. rows=" << M.Rows() << endl;
  //cerr << "no. cols=" << M.Cols() << endl;
  int i, j;
  for (i = 0; i < M.Rows (); i++) {
    for (j = 0; j < M.Cols (); j++) {
      file2 << M (i, j);
      file2 << endl;
    }
  }
  file2.close ();
}

void
write_matrix2 (irtkMatrix M)
{
  ofstream file2 ("matrixerror2", ios::out | ios::app);
  //cerr << "no. rows=" << M.Rows() << endl;
  //cerr << "no. cols=" << M.Cols() << endl;
  int i, j;
  for (i = 0; i < M.Rows (); i++) {
    for (j = 0; j < M.Cols (); j++) {
      file2 << M (i, j);
      file2 << endl;
    }
  }
  file2.close ();
}

int
main (int argc, char **argv)
{

  int i, j, k;
  int iNoOfLandmarks, iNoOfLandmarks_1, iNoOfLandmarks_2;
  int iNoOfDatasets, NoOfDatasetsAdded;
  double p[3];
  float scalarMean_1, scalarMean_2, Std_1, Std_2;
  float fTotalVar, fNorm;
  double nipalsTol = 0.00000000001;
  float norm;
  int nipalsIterations;

  irtkMatrix M, M_1, M_2;
  irtkMatrix T, T_1, T_2;
  irtkMatrix Eigenvector, Eigenvector_1, Eigenvector_2;
  irtkMatrix B_Eigenvector, data1_reduced, data2_reduced;
  irtkMatrix C_11, C_12, C_21, C_22;
  irtkMatrix U, V;
  irtkMatrix Finaleigenvector, Finaleigenvector_1, Finaleigenvector_2;
  irtkMatrix D_1, D_2, D;
  irtkMatrix M_1_predict, M_2_actual;
  irtkMatrix tempMatrix;

  irtkVector MeanShape, MeanShape_1, MeanShape_2, Marg_MeanShape_2;
  irtkVector Eigenvalues, Eigenvalues_1, Eigenvalues_2;

  scalarMean_1 = scalarMean_2 = Std_1 = Std_2 = 0.0;
  iNoOfLandmarks_1 = iNoOfLandmarks_2 = 0;

  if ((argc < 5)) {
    cout <<
      "usage: wpls_ibim [sub0_struc1] ... [subN-1_struc1] [sub0_struc2] ... [subN-1_struc2] [IndexOfSubjectLeftOut] [modesToUse] [maxiterations] [wplsurface.vtk] [meansurface.vtk] " << endl;
    // NumberOfSubjectLeftOut is indexed from zero
    // The other subjects are used to predict that subject's structre 2 from its structure 1.
    // Based on the data obtained from the other subjects.
    // Output file for the errors.
    return 1;
  }

  // total no of datasets, including set that will be predicted
  int TNoOfDatasets = (argc - 6) / 2;
  //no of datasets used for training
  iNoOfDatasets = TNoOfDatasets - 1;

  // get dataset number we will leave out, this is indexed from 0
  int leaveout = atoi (argv[2 * TNoOfDatasets + 1]);
  int noOfModes = atoi(argv[2 * TNoOfDatasets + 2]);
  int maxits= atoi (argv[2 * TNoOfDatasets + 3]);
  ///////////////////////////////////////////////////////////////////

  // Read all landmarks_1 of non-predicted data in matrix M_1, need to loop
  // over all datasets including predicted
  NoOfDatasetsAdded = 0;
  for (i = 0; i < TNoOfDatasets; i++) {

    data1    = read (argv[i + 1]);
    points1  = data1->GetPoints ();
    vectors1 = data1->GetPointData ()->GetVectors ();

    if (i == 0) {

      iNoOfLandmarks_1 = points1->GetNumberOfPoints ();

      M_1_predict.Initialize (3 * iNoOfLandmarks_1, 1);
      M_1.Initialize (3 * iNoOfLandmarks_1, iNoOfDatasets);

      MeanShape_1.Initialize (3 * iNoOfLandmarks_1);
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++) {
	MeanShape_1 (j) = 0.0;
      }

    } else if ( (points1->GetNumberOfPoints ()) != iNoOfLandmarks_1) {
      cerr << "Datasets must contain the same number of landmarks" << endl;
      exit (1);
    }

    if (i != leaveout){
      // include only non-predicted data
      for (j = 0; j < iNoOfLandmarks_1; j++) {
	points1->GetPoint (j, p);

	M_1 (3 * j    , NoOfDatasetsAdded) = p[0];
	M_1 (3 * j + 1, NoOfDatasetsAdded) = p[1];
	M_1 (3 * j + 2, NoOfDatasetsAdded) = p[2];

	scalarMean_1 += p[0];
	scalarMean_1 += p[1];
	scalarMean_1 += p[2];

	MeanShape_1 (3 * j    ) += p[0];
	MeanShape_1 (3 * j + 1) += p[1];
	MeanShape_1 (3 * j + 2) += p[2];
      }
      NoOfDatasetsAdded++;
    } else {
      cout << "LOO predictor: " << argv[i + 1] << endl;
      // read in data we will predict from
      for (j = 0; j < iNoOfLandmarks_1; j++) {
	points1->GetPoint (j, p);
	M_1_predict (3 * j    , 0) = p[0];
	M_1_predict (3 * j + 1, 0) = p[1];
	M_1_predict (3 * j + 2, 0) = p[2];
      }
    }
  } // Loop reading datasets 1.



  // Read all landmarks_2 in matrix M_2
  NoOfDatasetsAdded = 0;
  for (i = 0; i < TNoOfDatasets; i++) {

    data2    = read (argv[i + TNoOfDatasets + 1]);
    points2  = data2->GetPoints ();
    vectors2 = data2->GetPointData ()->GetVectors ();

    if (i == 0) {

      iNoOfLandmarks_2 = points2->GetNumberOfPoints ();

      M_2_actual.Initialize (3 * iNoOfLandmarks_2, 1);
      M_2.Initialize (3 * iNoOfLandmarks_2, iNoOfDatasets);

      MeanShape_2.Initialize (3 * iNoOfLandmarks_2);
      for (j = 0; j < 3 * iNoOfLandmarks_2; j++) {
	MeanShape_2 (j) = 0.0;
      }

    } else {
      if ((points2->GetNumberOfPoints ()) != iNoOfLandmarks_2) {
	cerr << "Datasets must contain the same number of landmarks" << endl;
	exit (1);
      }
    }

    if (i != leaveout){
      // include only non-predicted data
      for (j = 0; j < iNoOfLandmarks_2; j++) {
	points2->GetPoint (j, p);

	M_2 (3 * j    , NoOfDatasetsAdded) = p[0];
	M_2 (3 * j + 1, NoOfDatasetsAdded) = p[1];
	M_2 (3 * j + 2, NoOfDatasetsAdded) = p[2];

	scalarMean_2 += p[0];
	scalarMean_2 += p[1];
	scalarMean_2 += p[2];

	MeanShape_2 (3 * j) += p[0];
	MeanShape_2 (3 * j + 1) += p[1];
	MeanShape_2 (3 * j + 2) += p[2];
      }
      NoOfDatasetsAdded++;
    } else {
      cout << "LOO predicted: " << argv[i + TNoOfDatasets + 1] << endl;
      for (j = 0; j < iNoOfLandmarks_2; j++) {
	points2->GetPoint (j, p);

	M_2_actual (3 * j    , 0) = p[0];
	M_2_actual (3 * j + 1, 0) = p[1];
	M_2_actual (3 * j + 2, 0) = p[2];
      }
    }
  } // Loop reading datasets 2

  //////////////////////////////////////////////////////////////

  MeanShape_1  = MeanShape_1 * (1.0 / iNoOfDatasets);
  MeanShape_2  = MeanShape_2 * (1.0 / iNoOfDatasets);

  scalarMean_1 = scalarMean_1 / (iNoOfDatasets * iNoOfLandmarks_1 * 3);
  scalarMean_2 = scalarMean_2 / (iNoOfDatasets * iNoOfLandmarks_2 * 3);

  iNoOfLandmarks = iNoOfLandmarks_1 + iNoOfLandmarks_2;

  MeanShape.Initialize (3 * iNoOfLandmarks);

  for (i = 0; i < 3 * iNoOfLandmarks_1; i++){
    MeanShape (i) = MeanShape_1 (i);
  }

  for (i = 3 * iNoOfLandmarks_1; i < 3 * iNoOfLandmarks; i++){
    MeanShape (i) = MeanShape_2 (i - 3 * iNoOfLandmarks_1);
  }

  //Perform mean-centering. Do PCA on structures. Subtract the mean.
  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfLandmarks_1; j++) {
      M_1 (3 * j    , i) -= MeanShape_1 (3 * j    );
      M_1 (3 * j + 1, i) -= MeanShape_1 (3 * j + 1);
      M_1 (3 * j + 2, i) -= MeanShape_1 (3 * j + 2);
    }
  }

  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfLandmarks_2; j++) {
      M_2 (3 * j    , i) -= MeanShape_2 (3 * j    );
      M_2 (3 * j + 1, i) -= MeanShape_2 (3 * j + 1);
      M_2 (3 * j + 2, i) -= MeanShape_2 (3 * j + 2);
    }
  }

  ///////////////////// IS THIS USED???
  M.Initialize (3 * (iNoOfLandmarks_1 + iNoOfLandmarks_2), iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfLandmarks_1; j++) {
      M(3 * j    , i) = M_1 (3 * j    , i);
      M(3 * j + 1, i) = M_1 (3 * j + 1, i);
      M(3 * j + 2, i) = M_1 (3 * j + 2, i);
    }
  }

  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = iNoOfLandmarks_1; j < iNoOfLandmarks_1 + iNoOfLandmarks_2; j++) {
      M(3 * j    , i) = M_2 (3 * (j - iNoOfLandmarks_1)    , i);
      M(3 * j + 1, i) = M_2 (3 * (j - iNoOfLandmarks_1) + 1, i);
      M(3 * j + 2, i) = M_2 (3 * (j - iNoOfLandmarks_1) + 2, i);
    }
  }

  //cout << "allocated meanshape" << endl;

  // need to do PCA of M_1 and M_2 separately first

  // Form matrix T = (1/iNoOfDatasets) M^T * M which has dimensions
  // iNoOfDatasets x iNoOfDatasets
  T_1.Initialize (iNoOfDatasets, iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      T_1 (i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_1; k++){
	T_1 (i, j) += M_1 (k, i) * M_1 (k, j);
      }
    }
  }

  T_1 /= iNoOfDatasets;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea_1 (iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++) {
    for (int j = 0; j < iNoOfDatasets; j++) {
      ea_1.Matrix (i, j) = T_1 (i, j);
    }
  }
  ea_1.DecrSortEigenStuff ();

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M*M^T.
  // Done by M * evecs(T).
  tempMatrix.Initialize(iNoOfDatasets, iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++){
    for (j = 0; j < iNoOfDatasets; j++) {
      tempMatrix(i, j) = ea_1.Eigenvector(i, j);
    }
  }

  Eigenvector_1 = M_1 * tempMatrix;

  Eigenvalues_1.Initialize (iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++) {
    Eigenvalues_1 (i) = ea_1.Eigenvalue (i);
  }

  fTotalVar = 0;
  for (i = 0; i < iNoOfDatasets; i++) {
    fTotalVar += ea_1.Eigenvalue (i);
  }

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets; i++) {
    fNorm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks_1; j++) {
      fNorm += Eigenvector_1 (j, i) * Eigenvector_1 (j, i);
    }
    fNorm = sqrt (fNorm);

    if (100 * ea_1.Eigenvalue (i) / fTotalVar > MIN_NORM) { //if (1){
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++) {
	Eigenvector_1 (j, i) /= fNorm;
      }
    } else {
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++) {
	Eigenvector_1 (j, i) = 0;
      }
      Eigenvalues_1 (i) = 0;
    }
  }

  // End of 1st PCA.
  //////////////////////////////////////////////////////////

  T_2.Initialize (iNoOfDatasets, iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      T_2 (i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_2; k++)
	T_2 (i, j) += M_2 (k, i) * M_2 (k, j);
    }
  }
  T_2 /= iNoOfDatasets;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea_2 (iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      ea_2.Matrix (i, j) = T_2 (i, j);
    }
  }
  ea_2.DecrSortEigenStuff ();

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M*M^T.
  // Done by M * evecs(T).
  for (i = 0; i < iNoOfDatasets; i++){
    for (j = 0; j < iNoOfDatasets; j++) {
      tempMatrix(i, j) = ea_2.Eigenvector(i, j);
    }
  }
  Eigenvector_2 = M_2 * tempMatrix;

  Eigenvalues_2.Initialize (iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++) {
    Eigenvalues_2 (i) = ea_2.Eigenvalue (i);
  }

  fTotalVar = 0;

  for (i = 0; i < iNoOfDatasets; i++) {
    fTotalVar += ea_2.Eigenvalue (i);
  }

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets; i++) {
    fNorm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks_2; j++) {
      fNorm += Eigenvector_2 (j, i) * Eigenvector_2 (j, i);
    }
    fNorm = sqrt (fNorm);

    if (100 * ea_2.Eigenvalue (i) / fTotalVar > MIN_NORM) { //if (1){
      for (j = 0; j < 3 * iNoOfLandmarks_2; j++) {
	Eigenvector_2 (j, i) /= fNorm;
      }
    } else {
      for (j = 0; j < 3 * iNoOfLandmarks_2; j++) {
	Eigenvector_2 (j, i) = 0;
      }
      Eigenvalues_2 (i) = 0;
    }
  }

  // End of 2nd PCA
  /////////////////////////////////////

  //Now convert x-u data into pca coordinates using only first noOfModes eigenvectors

  // Get the reduced set of eigenvectors and find the components of the
  // data w.r.t. this set.
  tempMatrix = Eigenvector_1(0, 0, 3 * iNoOfLandmarks_1, noOfModes);
  tempMatrix.Transpose();
  data1_reduced = tempMatrix * M_1;

  tempMatrix = Eigenvector_2(0, 0, 3 * iNoOfLandmarks_2, noOfModes);
  tempMatrix.Transpose();
  data2_reduced= tempMatrix * M_2;

  // Means and sds calculated across all subjects and components to give a
  // single scalar value for each dataset. Req of NIPALS


  
  // Subtract the mean and calculate stds
  scalarMean_1 = scalarMean_2 = 0.0;
  for (i = 0; i < noOfModes; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      scalarMean_1 += data1_reduced (i, j);
      scalarMean_2 += data2_reduced (i, j);
    }
  }

  scalarMean_1 /= (noOfModes * iNoOfDatasets);
  scalarMean_2 /= (noOfModes * iNoOfDatasets);

  Std_1 = Std_2 = 0.0;
  for (i = 0; i < noOfModes; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      data1_reduced (i, j) -= scalarMean_1;
      Std_1 = Std_1 + data1_reduced (i, j) * data1_reduced (i, j);

      data2_reduced (i, j) -= scalarMean_2;
      Std_2 = Std_2 + data2_reduced (i, j) * data2_reduced (i, j);
    }
  }

  Std_1 = sqrt (Std_1 / (iNoOfDatasets * noOfModes));
  Std_2 = sqrt (Std_2 / (iNoOfDatasets * noOfModes));

  for (i = 0; i < noOfModes; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      data1_reduced (i, j) /= Std_1;
      data2_reduced (i, j) /= Std_2;
    }
  }
  cerr << "data1_reduced=" << endl;
  //print_matrix(data1_reduced);
  cerr << "data2_reduced=" << endl;
  //print_matrix(data2_reduced);
  
  // Now use NIPALS algorithm

  // Need to be in a form where each row is an observation (shape) to use
  // Nick's algorithm
  
  irtkMatrix X = data1_reduced;
  irtkMatrix Y = data2_reduced;
  X.Transpose ();
  Y.Transpose ();

  irtkMatrix A;
  irtkMatrix f;
  irtkMatrix A_T;
  irtkMatrix W_k_T;
  irtkMatrix P_k;
  irtkMatrix C;
  irtkMatrix Q_T;

  //irtkMatrix Q (noOfModes, iNoOfDatasets); old routine
  irtkMatrix Q (noOfModes, maxits);
  Q=Q*0.0;
  irtkMatrix Q_k (noOfModes, 1);

  irtkMatrix Q_k_T;
  irtkMatrix P_k_T;

  // X' * Y
  A = data1_reduced * Y;

  // X' * X
  M = data1_reduced * X;

  irtkMatrix B (noOfModes, noOfModes);
  B.Ident ();

  //irtkMatrix W (noOfModes, iNoOfDatasets); old routine
  irtkMatrix W (noOfModes, maxits);
  W=W*0.0;
  irtkMatrix W_k (noOfModes, 1);
  irtkMatrix A_TmultA;
  irtkMatrix firstevector(noOfModes,1);
  irtkMatrix eigcheck(noOfModes,1);
  // CODE TO ADD TO, REMEMBER ALSO QL ALGORITHM!!!!!!!!!
  float firstevalue;
  float ** choltest = matrix(1,noOfModes,1,noOfModes);
  float * choltest_vec = vector(1,noOfModes);
  float simcheck;
  // iterations chosen for nipals is maxits here - can choose other, depending on how many factors wanted.
  
  for (k = 0; k < maxits; k++){

    irtkEigenAnalysis ea_1 (noOfModes);
    A_T = A;
    A_T.Transpose ();

    A_TmultA = A_T * A;

    //CHECK SYMETTRY 
    simcheck = 0;
    for (i = 0; i < A_TmultA.Cols (); i++) {
      for (j = 0; j < A_TmultA.Cols (); j++) {
       simcheck+= fabs((A_TmultA)(i,j)-(A_TmultA)(j,i));
	 }
    }
    cout << "simcheck=" << simcheck << endl;
    //choldc(choltest,noOfModes,choltest_vec);


    for (i = 0; i < A.Cols (); i++) {
      for (j = 0; j < A.Cols (); j++) {
	//if (fabs((A_TmultA)(i,j))<1.0e-10) (A_TmultA)(i,j)=0.0;
	ea_1.Matrix (i, j) = (A_TmultA) (i, j);
      }
    }
    //print_matrix(A_TmultA);
    ea_1.DecrSortEigenStuff ();
    if (ea_1.error && 0x00000004 ){ //irtkEigenAnalysis::ql_exceeded
       cout << "Breaking out of nipals ... QL iterations exceeded." << endl;
       break;
     }

    //kth column of Q is last column of eigenvectors NEED TO INITIALIZE Q AND Q_k!
    for (i = 0; i < noOfModes; i++){
      Q (i, k) = Q_k (i, 0) = ea_1.Eigenvector (i, 0);
    }

    //CHECK EIGENVECTOR IS EIGENVECTOR! !!!!!!!!!!!!!!
    for (i = 0; i < noOfModes; i++) firstevector(i,0) = ea_1.Eigenvector (i, 0);
    firstevalue=ea_1.Eigenvalue(0);
    eigcheck= A_TmultA * firstevector - (firstevector*firstevalue);
    norm = 0.0;
    for (i = 0; i < noOfModes; i++) norm+=fabs(eigcheck(i,0));
    cout << "eigcheck= " << norm << endl;
    //kth column of W is B*A* kth column of Q NEED TO INITIALIZE W!
    
    norm = 0.0;
    tempMatrix = B * A * Q_k;
    for (i = 0; i < noOfModes; i++){
      W (i, k) = W_k (i, 0) = tempMatrix(i, 0);
    }

    for (i = 0; i < noOfModes; i++){
      norm += W (i, k) * W (i, k);
    }
    //cout << norm << endl;
    for (i = 0; i < noOfModes; i++){
      W (i, k) = W_k (i, 0) = W_k (i, 0) / sqrt (norm);
    }

    W_k_T = W_k;
    W_k_T.Transpose ();

    f = W_k_T * M * W_k;
    P_k = M * W_k / (f (0, 0));

    Q_k = A_T * W_k / (f (0, 0));;

    for (i = 0; i < noOfModes; i++){
      Q (i, k) = Q_k (i, 0);
    }

    Q_k_T = Q_k;
    Q_k_T.Transpose ();

    P_k_T = P_k;
    P_k_T.Transpose ();

    A = A - (P_k * Q_k_T * (f (0, 0)));
    M = M - (P_k * P_k_T * (f (0, 0)));

    B = B - W_k * P_k_T;

    if (sum_abs_matrix (M) < nipalsTol) 
    {
     cout << "M converged" << endl;
     break;
    }

    if (sum_abs_matrix (B) < nipalsTol)
    {
     cout << "B converged" << endl;
     break;
    }

    if (sum_abs_matrix (A) < nipalsTol)
    {
     cout << "A converged" << endl;
     break;
    }
    cout << "iteration=" << k << "matrix M norm=" << sum_abs_matrix (M)/(noOfModes*noOfModes) << "matrix B norm=" << sum_abs_matrix (B)/(noOfModes*noOfModes)
<< "matrix A norm=" <<sum_abs_matrix (A)/(noOfModes*noOfModes) <<"matrix AAT norm=" <<sum_abs_matrix (A_TmultA)/(noOfModes*noOfModes)
<< endl;
  }
  // End of main NIPALS loop.
 
  nipalsIterations = k;

  Q_T = Q;
  Q_T.Transpose ();

  C = W * Q_T;

  // 
  // Dataset of structure 1 is input. Reduce to its pca1 form. Then use C
  // and the mean1, std1 to give an output, then use mean2, std2 to give
  // this output in pca2 form. Then transform back to original basis. Could
  // in theory write C, the means, the stds and the eigenvectors into a
  // file and then just use this for the prediction.

  // Keep same as it is. Last argument, which is mean at the moment, can be
  // a number representing which pair from the datasets is to be the left
  // out one.

  // change M_1 and M_2 accordingly so that this dataset is not used to
  // create the model.  use C etc. as above to predict structure2 from
  // structure1 and then compare prediction with actual structure which was
  // read in originally!  perhaps use this format to do a cca testing as
  // well. note that there will be errors associated with reduction to pca,
  // as well as those produced by the wpls method itself

  irtkMatrix X_predict = M_1_predict;

  // firstly put into pca basis
  for (i = 0; i < 3 * iNoOfLandmarks_1; i++){
    //subtract Mean
    X_predict (i, 0) = X_predict (i, 0) - MeanShape_1 (i);
  }

  //Convert to pca coordinates using only first noOfModes eigenvectors
  irtkMatrix X_predict_pca;
  tempMatrix = Eigenvector_1(0, 0, 3 * iNoOfLandmarks_1, noOfModes);
  tempMatrix.Transpose();
  X_predict_pca = tempMatrix * X_predict;

  // calculate error due to pca modelling for structure 1
  irtkMatrix X_predict_pca_back (3 * iNoOfLandmarks_1, 1);

  for (i = 0; i < 3 * iNoOfLandmarks_1; i++) {
    X_predict_pca_back (i, 0) = 0.0;
    for (k = 0; k < noOfModes; k++){
      X_predict_pca_back (i, 0) += Eigenvector_1 (i, k) * X_predict_pca (k, 0);
    }
  }

  // dont add mean back, as we compare with X_predict which also had mean subtracted
  float pcaX_error = 0;
  for (i = 0; i < 3 * iNoOfLandmarks_1; i++){
    pcaX_error += (X_predict (i, 0) - X_predict_pca_back (i, 0)) *
                  (X_predict (i, 0) - X_predict_pca_back (i, 0));
  }

  pcaX_error = pcaX_error / (iNoOfLandmarks_1);

  // now subtract matrix mean and divide by std
  X_predict_pca -= (scalarMean_1);
  X_predict_pca /= (Std_1);

  X_predict_pca.Transpose ();
  // input structure is now in right form for prediction
  //cout << "predicting Y" << endl;

  // HERE IS THE ACTUAL PREDICTION!

  irtkMatrix Y_predicted_pca = X_predict_pca * C;
  //cout << "predicted Y" << endl;
  // now convert this back
  Y_predicted_pca *= (Std_2);
  Y_predicted_pca += (scalarMean_2);

  Y_predicted_pca.Transpose ();

  //convert actual Y into pca coordinates for comparison with ypredictedpca
  //SO THAT WE CAN COMPARE THE DIFFERENCE WITH THE PREDICTED Y IN PCA
  //COORDINATES.

  irtkMatrix M_2_actual_nomean (3 * iNoOfLandmarks_2, 1);

  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    M_2_actual_nomean (i, 0) = M_2_actual (i, 0) - MeanShape_2 (i);
  }

  irtkMatrix Y_actual_pca (noOfModes, 1);

  for (i = 0; i < noOfModes; i++) {
    Y_actual_pca (i, 0) = 0.0;

    for (k = 0; k < 3 * iNoOfLandmarks_2; k++){
      Y_actual_pca (i, 0) += Eigenvector_2 (k, i) * M_2_actual_nomean (k, 0);
    }
  }

//   cout << "struct2pca=" << endl;
//   print_matrix (Y_actual_pca);

  // convert back into normal coordinates, to calculate error of pca model, this gives a sort of lower bound for our predicted error
  irtkMatrix Y_actual_pca_back (3 * iNoOfLandmarks_2, 1);
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++) {
    Y_actual_pca_back (i, 0) = 0.0;

    for (k = 0; k < noOfModes; k++){
      Y_actual_pca_back (i, 0) += Eigenvector_2 (i, k) * Y_actual_pca (k, 0);
    }
  }

  // dont add mean back, as we compare with M_2_actual_nomean which also had mean subtracted
  float pcaY_error = 0;

  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    pcaY_error += (Y_actual_pca_back (i, 0) - M_2_actual_nomean (i, 0)) * 
                  (Y_actual_pca_back (i, 0) - M_2_actual_nomean (i, 0));
  }

  pcaY_error = pcaY_error / (iNoOfLandmarks_2);


  // prediction is now in pca form of structure 2. keep this because we might use it for evaluation
  irtkMatrix Y_predicted (3 * iNoOfLandmarks_2, 1);

  for (i = 0; i < 3 * iNoOfLandmarks_2; i++) {
    Y_predicted (i, 0) = 0.0;
    for (k = 0; k < noOfModes; k++){
      Y_predicted (i, 0) += Eigenvector_2 (i, k) * Y_predicted_pca (k, 0);
    }
  }

  // addback Mean
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    Y_predicted (i, 0) = Y_predicted (i, 0) + MeanShape_2 (i);
  }




  // now finally do comparison with actual shape! IN THE ORIGINAL SPACE ....
  float wplsError = 0.0;
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    wplsError += (Y_predicted (i, 0) - M_2_actual (i, 0)) * 
                 (Y_predicted (i, 0) - M_2_actual (i, 0));
  }

  wplsError = wplsError / (iNoOfLandmarks_2);

  // calc error of using mean of 2ndshape SQUARED DIFFERENCE BETWEEN MEAN 2 AND ACTUAL 2
  float meanError = 0.0;
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    meanError += (MeanShape_2 (i) - M_2_actual (i, 0)) *
                 (MeanShape_2 (i) - M_2_actual (i, 0));
  }

  meanError = meanError / (iNoOfLandmarks_2);

  cout << iNoOfDatasets << ",";
  cout << leaveout << ",";
  cout << iNoOfLandmarks_1 << ",";
  cout << iNoOfLandmarks_2 << ",";
  cout << nipalsIterations << ",";
  cout << wplsError << ",";
  cout << meanError << ",";
  cout << pcaX_error << ",";
  cout << pcaY_error << endl;
// Read surfaces


 vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
 surface_reader->SetFileName(argv[2*TNoOfDatasets]);
 surface_reader->Modified();
 surface_reader->Update();
 vtkPolyData *surface = surface_reader->GetOutput();


 for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    meanError += (MeanShape_2 (i) - M_2_actual (i, 0)) *
                 (MeanShape_2 (i) - M_2_actual (i, 0));
  }

 vtkFloatArray * scalars_meanerror = vtkFloatArray::New();
 vtkFloatArray * scalars_wplserror = vtkFloatArray::New();

 double coord[3];
 for (i = 0; i < surface->GetNumberOfPoints(); i++)
 {
   coord[0]=Y_predicted(3*i,0);
   coord[1]=Y_predicted(3*i+1,0);
   coord[2]=Y_predicted(3*i+2,0);
  surface->GetPoints()->SetPoint(i, coord);
 }
 surface->Modified();
 float value;
 for (i = 0; i < iNoOfLandmarks_2; i++){
    value= (Y_predicted(3*i,0) - M_2_actual (3*i, 0))*(Y_predicted(3*i,0) - M_2_actual (3*i, 0))+
    (Y_predicted(3*i+1,0) - M_2_actual (3*i+1, 0))*(Y_predicted(3*i+1,0) - M_2_actual (3*i+1, 0))+
    (Y_predicted(3*i+2,0) - M_2_actual (3*i+2, 0))*(Y_predicted(3*i+2,0) - M_2_actual (3*i+2, 0));
    value=sqrt(value);
    scalars_wplserror->InsertTuple1(i, value);
 }
 surface->GetPointData()->AddArray(scalars_wplserror);
 vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
 writer->SetInput(surface);
 writer->SetFileName(argv[2*TNoOfDatasets+4]);
 writer->Write();

 for (i = 0; i < surface->GetNumberOfPoints(); i++)
 {
   coord[0]=MeanShape_2(3*i);
   coord[1]=MeanShape_2(3*i+1);
   coord[2]=MeanShape_2(3*i+2);
  surface->GetPoints()->SetPoint(i, coord);
 }
 surface->Modified();
 for (i = 0; i < iNoOfLandmarks_2; i++){
    value= (MeanShape_2 (3*i) - M_2_actual (3*i, 0))*(MeanShape_2 (3*i) - M_2_actual (3*i, 0))
    +(MeanShape_2 (3*i+1) - M_2_actual (3*i+1, 0))*(MeanShape_2 (3*i+1) - M_2_actual (3*i+1, 0))+
    (MeanShape_2 (3*i+2) - M_2_actual (3*i+2, 0))*(MeanShape_2 (3*i+2) - M_2_actual (3*i+2, 0));
    value=sqrt(value);
    scalars_meanerror->InsertTuple1(i, value);
 }
 surface->GetPointData()->AddArray(scalars_meanerror);

 vtkPolyDataWriter *writer2 = vtkPolyDataWriter::New();
 writer2->SetInput(surface);
 writer2->SetFileName(argv[2*TNoOfDatasets+5]);
 writer2->Write();





}

#else

#include <irtkImage.h>

int
main (int argc, char *argv[])
{
  cout << argv[0] << " needs to be compiled with the contrib and VTK library"
    << endl;
}
#endif





  /*
     // calc error of using marginal mean of 2nd shape ie. mean of 2nd shape given test shape 1
     Marg_MeanShape_2.Initialize(3*iNoOfLandmarks_2);

     //marginal mean=MeanShape_2(i)+C_21*(C_11)^(-1)*(M_1_actual-MeanShape_1(i))

     C_11.Initialize(modeskept,modeskept);
     C_22.Initialize(modeskept,modeskept);
     C_12.Initialize(modeskept,modeskept);
     C_21.Initialize(modeskept,modeskept);

     // calculate covariance matrices for reduced_data
     C_11.Initialize(modeskept,modeskept);
     for (i = 0; i < modeskept; i++){
     for (j = 0; j < modeskept; j++){
     C_11(i, j) = 0.0;
     for (k = 0; k < iNoOfDatasets; k++) C_11(i, j) += data1_reduced_copy(i,k)*data1_reduced_copy(j,k);
     }
     }
     C_11 /= iNoOfDatasets;

     C_22.Initialize(modeskept,modeskept);
     for (i = 0; i < modeskept; i++){
     for (j = 0; j < modeskept; j++){
     C_22(i, j) = 0.0;
     for (k = 0; k < iNoOfDatasets; k++) C_22(i, j) += data2_reduced_copy(i,k)*data2_reduced_copy(j,k);
     }
     }
     C_22 /= iNoOfDatasets;

     C_12.Initialize(modeskept,modeskept);
     for (i = 0; i < modeskept; i++){
     for (j = 0; j < modeskept; j++){
     C_12(i, j) = 0.0;
     for (k = 0; k < iNoOfDatasets; k++) C_12(i, j) += data1_reduced_copy(i,k)*data2_reduced_copy(j,k);
     }
     }
     C_12 /= iNoOfDatasets;

     C_21.Initialize(modeskept,modeskept);
     for (i = 0; i < modeskept; i++){
     for (j = 0; j < modeskept; j++){
     C_21(i, j) = C_12(j, i);
     }
     }

     irtkMatrix X_predict_pca_copy(modeskept,1);
     //Convert input shape into x-u pca coordinates using only first modeskept eigenvectors
     for (i = 0; i <  modeskept; i++){
     X_predict_pca_copy(i,0)=0.0;
     for (k=0;k < 3*iNoOfLandmarks_1;k++) X_predict_pca_copy(i,0)+=Eigenvector_1(k,i)*X_predict(k,0);
     }


     irtkMatrix Z;
     C_11.Invert();
     //Z=C_21*C_11*X_predict_pca; Z now has marginal mean in pca x-u coords
     //X_predict_pca_copy.Transpose();
     Z=C_21*C_11*X_predict_pca_copy; // Z now has marginal mean in pca x-u coords
     // must convert this back to original coords

     for (i = 0; i < 3*iNoOfLandmarks_2; i++){
     Marg_MeanShape_2(i)=0.0;
     //for (k=0;k<modeskept;k++) Finaleigenvector_1(i,j)+=Eigenvector_1(i,k)*can_1.Eigenvector(k,j);
     for (k=0;k<modeskept;k++) Marg_MeanShape_2(i)+=Eigenvector_2(i,k)*Z(k,0);
     }
     for (i = 0; i < 3*iNoOfLandmarks_2; i++) Marg_MeanShape_2(i)+=MeanShape_2(i);
     for (i=0;i<3*iNoOfLandmarks_2;i++) error3=error3+(M_2_actual(i,0)-Marg_MeanShape_2(i))*(M_2_actual(i,0)-Marg_MeanShape_2(i));
     error3=error3/(iNoOfLandmarks_2);
   */

  //cout << error << " " << meanError << " " << error3 << " " << pcaX_error << " " << pcaY_error << endl;
  //ofstream file1 ("/vol/vipdata/users/ar17/brainproject/93_subjects/wpls/errors", ios::out | ios::app );
//   ofstream file1 (argv[2 * TNoOfDatasets + 2], ios::out | ios::app);
//   //file1 << leaveout << " " << error << " " << meanError << " " << error3 << " " << pcaX_error << " " << pcaY_error << " " << pcaY_error-error<< "\n";
//   file1 << leaveout << " " << error << " " << meanError << " " << pcaX_error <<
//     " " << pcaY_error << " " << pcaY_error - error << "\n";
//   file1.close ();

