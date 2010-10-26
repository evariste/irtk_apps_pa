#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <fstream.h>
#include <stdio.h>
#include <string.h>
#include <irtkImage.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>
#include <nr.h>
#include <nrutil.h>

// Minimum norm of an eigenvector
#define MIN_NORM 0.01

// vtk includes
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkPolyDataAlgorithm.h>

#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

vtkPointSet  *data1;
vtkPoints    *points1;
vtkDataArray *vectors1;

vtkPointSet  *data2;
vtkPoints    *points2;
vtkDataArray *vectors2;

char *dumbMeanOutputFile             = NULL;
char *marginalMeanOutputFile         = NULL;
char *inputDependentDataTemplateFile = NULL;
char *predictingDataOutputMatrix     = NULL;
char *dependentDataOutputMatrix      = NULL;


// Reads points / polys from a file.
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

  for (i = 0; i < 3; i++){
    is.getline(buffer, 256);
  }
  is >> buffer >> buffer;

  if (strcmp(buffer, "POLYDATA") == 0) {

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


// Writes points / polys to a file.
void write(char *file, vtkPointSet *pset)
{
  if (pset->IsA("vtkStructuredGrid")){
    vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
    writer->SetInput((vtkStructuredGrid *)pset);
    writer->SetFileName(file);
    writer->SetFileTypeToBinary();
    writer->Update();
  } else {
    if (pset->IsA("vtkUnstructuredGrid")){
      vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
      writer->SetInput((vtkUnstructuredGrid *)pset);
      writer->SetFileName(file);
      writer->SetFileTypeToBinary();
      writer->Update();
    } else {
      if (pset->IsA("vtkPolyData")){
        vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
        writer->SetInput((vtkPolyData *)pset);
        writer->SetFileName(file);
        writer->SetFileTypeToASCII();
        writer->Write();
      } else {
        cerr << "Unknown VTK data type" << endl;
        exit(1);
      }
    }
  }
}

//
void usage(){
  cout << "usage: marginalmean_ibim  N indexLeftOut pred_0 ... pred_n-1 dep_0 ... dep_n-1";
  cout << " dumbMeanOutputFile marginalMeanOutputFile <-options>" << endl;
  cout << "" << endl;
  cout << "N is the number of data sets in each group (predicting and dependent)" << endl;
  cout << "indexLeftOut refers to the subject to be left out, indexed from zero." << endl;

  cout << "Writes out dumb mean and marginal mean n for the dependent structure" << endl;
  cout << "of the left out subject." << endl;
  cout << "Joint distribution based on all data sets except left out one." << endl;

  cout << "<-modes> number     : How many modes to use when representing the data in" << endl;
  cout << "                      the PCA basis. Default = 13." << endl;

  cout << "<-matrix_pred> file : write the predicting data matrix to a file." << endl;
  cout << "<-matrix_dep> file  : write the dependent data matrix to a file." << endl;

  cout << "" << endl;
  exit(1);
}

//
int main(int argc, char **argv)
{
  int i, j, k;
  bool ok;
  int iNoOfLandmarks, iNoOfLandmarks_1, iNoOfLandmarks_2;
  int iNoOfDatasets;

  irtkMatrix M, M_1, M_2;
  irtkMatrix T, T_1, T_2;
  irtkMatrix Eigenvector, Eigenvector_1, Eigenvector_2;
  irtkMatrix data1_reduced, data2_reduced; 
  irtkMatrix C_11, C_12, C_21, C_22;
  irtkMatrix M_1_predict, M_2_actual;

  irtkVector MeanShape, MeanShape_1, MeanShape_2;
  irtkVector marginalMeanShape2;
  irtkVector Eigenvalues, Eigenvalues_1, Eigenvalues_2;

  double p[3];
  double coord[3];
  float error1, error2;

  // Number of modes to use, default value.
  int R = 13;

  if ((argc < 5)){
    usage();
  }

  int NoOfDatasetsAdded = 0;
  // Total no of datasets, including set that will be predicted.
  //   int TNoOfDatasets = (argc - 4) / 2;
  int TNoOfDatasets = atoi(argv[1]);
  argv++;
  argc--;

  // Get dataset number we will leave out, this is indexed from 0
  int leaveout = atoi( argv[1] );
  argv++;
  argc--;

  // No. of datasets used for training.
  iNoOfDatasets = TNoOfDatasets - 1;
  cerr << " There are " << iNoOfDatasets << " data sets for training" << endl;
  cerr << " The index left out is " << leaveout << endl;

  iNoOfLandmarks_1 = iNoOfLandmarks_2 = 0;

  // Read all landmarks_1 of non-predicted data in matrix M_1, need to loop
  // over all datasets including predicted
  for (i = 0; i < TNoOfDatasets; i++){

    if (i == leaveout){
      cout << "Leaving out " << argv[1] << endl;
    }

    data1    = read(argv[1]);
    argv++;
    argc--;

    points1  = data1->GetPoints();
    vectors1 = data1->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks_1 = points1->GetNumberOfPoints();

      M_1_predict.Initialize(3 * iNoOfLandmarks_1, 1);
      M_1.Initialize(3 * iNoOfLandmarks_1, iNoOfDatasets);
      MeanShape_1.Initialize(3 * iNoOfLandmarks_1);

      for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
        MeanShape_1(j) = 0.0;
      }

      cout << " There are " << iNoOfDatasets << " datasets with ";
      cout << points1->GetNumberOfPoints()   << " landmarks." << endl;
    } else if (points1->GetNumberOfPoints() != iNoOfLandmarks_1){
      cerr << "Datasets must contain the same number of landmarks" << endl;
      exit(1);
    }

    if (i != leaveout){

      // include only non-predicted data, i.e. the training data.
      for (j = 0; j < iNoOfLandmarks_1; j++){

        points1->GetPoint(j, p);

        M_1(3 * j    , NoOfDatasetsAdded) = p[0];
        M_1(3 * j + 1, NoOfDatasetsAdded) = p[1];
        M_1(3 * j + 2, NoOfDatasetsAdded) = p[2];

        // Mean shape, however is needed for PCA.
        MeanShape_1(3 * j    ) += p[0];
        MeanShape_1(3 * j + 1) += p[1];
        MeanShape_1(3 * j + 2) += p[2];
      }
      NoOfDatasetsAdded++;
    } else {

      // Read in data we will predict from, left-one-out subject.
      for (j = 0; j < iNoOfLandmarks_1; j++){
        points1->GetPoint(j, p);
        M_1_predict(3 * j    , 0) = p[0];
        M_1_predict(3 * j + 1, 0) = p[1];
        M_1_predict(3 * j + 2, 0) = p[2];
      }
    }
  }

  NoOfDatasetsAdded = 0;
  cerr << "Read predicting data:" << endl;
  // Read all landmarks_2 into matrix M_2.

  for (i = 0; i < TNoOfDatasets; i++){

    if (i == leaveout){
      cout << "Leaving out " << argv[1] << endl;
    }

    if (i == TNoOfDatasets - 1){
      // Use the last of the dependent data sets as a template for writing
      // the prediction, could be any of the dependents, the main thing is
      // to get the right connectivity.
      inputDependentDataTemplateFile = argv[1];
    }

    data2    = read(argv[1]);
    argv++;
    argc--;

    points2  = data2->GetPoints();
    vectors2 = data2->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks_2 = points2->GetNumberOfPoints();

      M_2_actual.Initialize(3 * iNoOfLandmarks_2, 1);
      M_2.Initialize(3 * iNoOfLandmarks_2, iNoOfDatasets);
      MeanShape_2.Initialize(3 * iNoOfLandmarks_2);

      for (j = 0; j < 3 * iNoOfLandmarks_2; j++){
        MeanShape_2(j) = 0.0;
      }

      cout << " There are " << iNoOfDatasets << " datasets with ";
      cout << iNoOfLandmarks_2 << " landmarks." << endl;
    } else {
      if ((points2->GetNumberOfPoints()) != iNoOfLandmarks_2){
        cerr << "Datasets must contain the same number of landmarks" << endl;
        exit(1);
      }
    }

    if (i != leaveout) {
      // include only non-predicted (training) data
      for (j = 0; j < iNoOfLandmarks_2; j++){

        points2->GetPoint(j, p);

        M_2(3 * j    , NoOfDatasetsAdded) = p[0];
        M_2(3 * j + 1, NoOfDatasetsAdded) = p[1];
        M_2(3 * j + 2, NoOfDatasetsAdded) = p[2];

        MeanShape_2(3 * j)   += p[0];
        MeanShape_2(3 * j + 1) += p[1];
        MeanShape_2(3 * j + 2) += p[2];
      }
      NoOfDatasetsAdded++;
    } else {
      for (j = 0; j < iNoOfLandmarks_2; j++){

        points2->GetPoint(j, p);

        M_2_actual(3 * j,     0) = p[0];
        M_2_actual(3 * j + 1, 0) = p[1];
        M_2_actual(3 * j + 2, 0) = p[2];
      }
    }
  }

  cerr << "read all landmarks" << endl;

  dumbMeanOutputFile = argv[1];
  argv++;
  argc--;
  marginalMeanOutputFile = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-matrix_pred") == 0)){
      argc--;
      argv++;
      predictingDataOutputMatrix = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-matrix_dep") == 0)){
      argc--;
      argv++;
      dependentDataOutputMatrix = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-modes") == 0)){
      argc--;
      argv++;
      R = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  cout << "Retaining " << R << " modes after PCA of the data." << endl;

  // Write debug data if needed.
  if (predictingDataOutputMatrix != NULL){
    M_1.Write(predictingDataOutputMatrix);
  }

  if (dependentDataOutputMatrix != NULL){
    M_2.Write(dependentDataOutputMatrix);
  }

  // Calculate the means.
  MeanShape_1 = MeanShape_1 * (1.0 / iNoOfDatasets);
  MeanShape_2 = MeanShape_2 * (1.0 / iNoOfDatasets);

  iNoOfLandmarks = iNoOfLandmarks_1 + iNoOfLandmarks_2;

  // Perform mean-centering and do PCA on structures.

  // Subtract the means.
  for (i = 0; i < iNoOfDatasets; i++){
    for (j = 0; j < iNoOfLandmarks_1; j++){
      M_1(3 * j,     i)  -= MeanShape_1(3 * j    );
      M_1(3 * j + 1, i)  -= MeanShape_1(3 * j + 1);
      M_1(3 * j + 2, i)  -= MeanShape_1(3 * j + 2);
    }
  }

  for (i = 0; i < iNoOfDatasets; i++){
    for (j = 0; j < iNoOfLandmarks_2; j++){
      M_2(3 * j,     i)  -= MeanShape_2(3 * j    );
      M_2(3 * j + 1, i)  -= MeanShape_2(3 * j + 1);
      M_2(3 * j + 2, i)  -= MeanShape_2(3 * j + 2);
    }
  }

  //
  // Need to do PCA of M_1 and M_2 separately first
  //

  // Form the 'small' covariance matrix 
  //   T = (1/iNoOfDatasets) M^T  *  M
  // which has dimensions
  //   iNoOfDatasets x iNoOfDatasets
  // This swaps the roles of variables and observations.

  T_1.Initialize(iNoOfDatasets, iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++){
    for (j = 0; j < iNoOfDatasets; j++){
      T_1(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_1; k++)
        T_1(i, j) += M_1(k, i) * M_1(k, j);
    }
  }

  T_1 /= iNoOfDatasets;

  // Compute the first eigen-decomposition.
  irtkEigenAnalysis ea_1(iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++){
    for(j = 0; j < iNoOfDatasets; j++){
      ea_1.Matrix(i, j) = T_1(i, j);
    }
  }
  ea_1.DecrSortEigenStuff();

  // Convert the eigenvectors of T to true eigenvectors of the covariance
  // matrix of the original data, the 'large' matrix:
  //   C = M * M^T.
  // This is due to ... if x is an eigenvector of T with eigenvalue c:
  //                T x = c x
  //      =>    M^T M x = c x
  //      =>  M M^T M x = c M x
  //      =>    M M^T y = c y
  //
  //      where y = M x
  //
  //      i.e. Mx is an eigenvector of the full covariance matrix of the
  //      original data (M * M^T)

  Eigenvector_1.Initialize(3 * iNoOfLandmarks_1, iNoOfDatasets);
  Eigenvector_1 *= 0.0;

  for (i = 0; i < iNoOfDatasets ; i++){
    for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
      for (k = 0; k < iNoOfDatasets; k++){
        Eigenvector_1(j, i) += M_1(j, k) * ea_1.Eigenvector(k, i);
      }
    }
  }

  Eigenvalues_1.Initialize(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets ; i++){
    Eigenvalues_1(i) = ea_1.Eigenvalue(i);
  }

  float fTotalVar   = 0;
  float fCummulated = 0;

  for (i = 0; i < iNoOfDatasets; i++){
    fTotalVar += ea_1.Eigenvalue(i);
  }

  // Some reporting.
  for (i = 0; i < R; i++){
    fCummulated += ea_1.Eigenvalue(i);
  }

  cout << "The first " << R << " modes of the PCA explain " << 100.0 * fCummulated / fTotalVar;
  cout << "% of the total variance of the predicting data." << endl;

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets ; i++){
    float fNorm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
      fNorm += Eigenvector_1(j, i) * Eigenvector_1(j, i);
    }

    fNorm = sqrt(fNorm);

    if (100 * ea_1.Eigenvalue(i) / fTotalVar > MIN_NORM){
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
        Eigenvector_1(j, i) /= fNorm;
      }
    } else {
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
        Eigenvector_1(j, i) = 0;
      }
      Eigenvalues_1(i) = 0;
    }
  }

  cerr << "done 1st PCA" << endl;

  //
  // Do PCA on second set of shapes, a bit of repetition ...
  //

  T_2.Initialize(iNoOfDatasets, iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; i++){
    for (j = 0; j < iNoOfDatasets; j++){
      T_2(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_2; k++)
        T_2(i, j) += M_2(k, i) * M_2(k, j);
    }
  }
  T_2 /= iNoOfDatasets;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea_2(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++){
    for(j = 0; j < iNoOfDatasets; j++){
      ea_2.Matrix(i, j) = T_2(i, j);
    }
  }
  ea_2.DecrSortEigenStuff();

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M * M^T.
  Eigenvector_2.Initialize(3 * iNoOfLandmarks_2, iNoOfDatasets);
  Eigenvector_2 *= 0.0;
  for (i = 0; i < iNoOfDatasets ; i++){
    for (j = 0; j < 3 * iNoOfLandmarks_2; j++){
      for (k = 0; k < iNoOfDatasets; k++){
        Eigenvector_2(j, i) += M_2(j, k) * ea_2.Eigenvector(k, i);
      }
    }
  }

  Eigenvalues_2.Initialize(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets ; i++){
    Eigenvalues_2(i) = ea_2.Eigenvalue(i);
  }

  fTotalVar   = 0;
  fCummulated = 0;

  for (i = 0; i < iNoOfDatasets; i++){
    fTotalVar += ea_2.Eigenvalue(i);
  }

  // Some reporting.
  for (i = 0; i < R; i++){
    fCummulated += ea_2.Eigenvalue(i);
  }

  cout << "The first " << R << " modes of the PCA explain " << 100.0 * fCummulated / fTotalVar;
  cout << "% of the total variance of the dependent data." << endl;

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets ; i++){
    float fNorm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks_2; j++){
      fNorm += Eigenvector_2(j, i) * Eigenvector_2(j, i);
    }

    fNorm = sqrt(fNorm);

    if (100 * ea_2.Eigenvalue(i) / fTotalVar > MIN_NORM){
      for (j = 0; j < 3 * iNoOfLandmarks_2; j++){
        Eigenvector_2(j, i) /= fNorm;
      }
    } else {
      for (j = 0; j < 3 * iNoOfLandmarks_2; j++){
        Eigenvector_2(j, i) = 0;
      }
      Eigenvalues_2(i) = 0;
    }
  }

  cerr << "Done 2nd PCA" << endl;

  //
  // Now convert x-u data into PCA coordinates using only first R
  // eigenvectors (Default value is 13)
  //

  if (R > iNoOfDatasets){
    cerr << "Number of modes required exceeds number available." << endl;
    exit(1);
  }

  data1_reduced.Initialize(R, iNoOfDatasets);
  for (i = 0; i <  R; i++){
    for (j = 0; j < iNoOfDatasets; j++){
      data1_reduced(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_1; k++){
        data1_reduced(i, j) += Eigenvector_1(k, i) * M_1(k, j);
      }
    }
  }

  data2_reduced.Initialize(R, iNoOfDatasets);
  for (i = 0; i <  R; i++){
    for (j = 0; j < iNoOfDatasets; j++){
      data2_reduced(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks_2; k++){
        data2_reduced(i, j) += Eigenvector_2(k, i) * M_2(k, j);
      }
    }
  }

  // Calculate error of using marginal mean of 2nd shape ie. mean of 2nd
  // shape given test shape 1.
  marginalMeanShape2.Initialize(3 * iNoOfLandmarks_2);

  //
  // HERE'S THE RUB!
  //

  // marginal mean = MeanShape_2(i) + C_21 * (C_11)^(-1) * (M_1_actual - MeanShape_1(i))

  // In the above MeanShape_2 and MeanShape_1 should be zero if we are
  // working in the pca bases.  Which explains why they do not actually
  // appear in ** below

  // Everything here is done in the pca bases.

  C_11.Initialize(R, R);
  C_22.Initialize(R, R);
  C_12.Initialize(R, R);
  C_21.Initialize(R, R);

  // Calculate covariance matrices for reduced_data
  C_11.Initialize(R, R);
  for (i = 0; i < R; i++){
    for (j = 0; j < R; j++){
      C_11(i, j) = 0.0;
      for (k = 0; k < iNoOfDatasets; k++){
        C_11(i, j) += data1_reduced(i, k) * data1_reduced(j, k);
      }
    }
  }
  C_11 /= iNoOfDatasets;

  C_22.Initialize(R, R);
  for (i = 0; i < R; i++){
    for (j = 0; j < R; j++){
      C_22(i, j) = 0.0;
      for (k = 0; k < iNoOfDatasets; k++){
        C_22(i, j) += data2_reduced(i, k) * data2_reduced(j, k);
      }
    }
  }
  C_22 /= iNoOfDatasets;
  // C_22 not actually used after this point.

  C_12.Initialize(R, R);
  for (i = 0; i < R; i++){
    for (j = 0; j < R; j++){
      C_12(i, j) = 0.0;
      for (k = 0; k < iNoOfDatasets; k++){
        C_12(i, j) += data1_reduced(i, k) * data2_reduced(j, k);
      }
    }
  }
  C_12 /= iNoOfDatasets;

  C_21.Initialize(R, R);
  for (i = 0; i < R; i++){
    for (j = 0; j < R; j++){
      C_21(i, j) = C_12(j, i);
    }
  }

  // Copy the predictor variables for the left out subject.
  irtkMatrix X_predict = M_1_predict;

  // Express the predictor for the left out subject in the predictor data PCA basis.

  //subtract Mean
  for (i = 0; i < 3 * iNoOfLandmarks_1; i++){

    X_predict(i, 0) = X_predict(i, 0) - MeanShape_1(i);
  }

  //Convert to PCA coordinates using first R eigenvectors.
  irtkMatrix X_predict_pca(R, 1);

  for (i = 0; i <  R; i++){
    X_predict_pca(i, 0) = 0.0;
    for (k = 0; k < 3 * iNoOfLandmarks_1; k++){
      X_predict_pca(i, 0) += Eigenvector_1(k, i) * X_predict(k, 0);
    }
  }

  // Now make the prediction - The result will be in the PCA basis of the
  // dependent data.
  irtkMatrix Z;

  //
  // HERE IS ** referred to earlier.
  //
  C_11.Invert();
  Z = C_21 * C_11 * X_predict_pca;
  //  RxR  * RxR  * R*1  <--- Sizes


  // Z now contains the marginal mean with respect to the PCA basis of the
  // dependent data.  Must convert this back to original coords

  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    marginalMeanShape2(i) = 0.0;
    for (k = 0; k < R; k++){
      marginalMeanShape2(i) += Eigenvector_2(i, k) * Z(k, 0);
    }
  }

  // Add the mean back in.
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    marginalMeanShape2(i) += MeanShape_2(i);
  }

  //
  // Some error checking.
  //

  error1 = 0.0;
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    error1 = error1 + (M_2_actual(i, 0) - marginalMeanShape2(i)) * (M_2_actual(i, 0) - marginalMeanShape2(i));
  }
  error1 /= ((double) iNoOfLandmarks_2);

  error2 = 0.0;
  for (i = 0; i < 3 * iNoOfLandmarks_2; i++){
    error2 = error2 + (M_2_actual(i, 0) - MeanShape_2(i)) * (M_2_actual(i, 0) - MeanShape_2(i));
  }
  error2 /= ((double) iNoOfLandmarks_2);

  cerr << "error of marginal mean = " << error1 << endl;
  cerr << "error of normal mean   = " << error2 << endl;

  //
  // write surface to file.
  //

  // File read here is the last of the dependent data sets, could be any of
  // the dependents, the main thing is to get the right connectivity.
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(inputDependentDataTemplateFile);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  for (i = 0; i < surface->GetNumberOfPoints(); i++)
  {
    coord[0] = marginalMeanShape2(3 * i    );
    coord[1] = marginalMeanShape2(3 * i + 1);
    coord[2] = marginalMeanShape2(3 * i + 2);
    surface->GetPoints()->SetPoint(i, coord);
  }
  surface->Modified();

  //Update the normals to reflect the new points.
  vtkPolyDataNormals *marginalNormals = vtkPolyDataNormals::New();
  marginalNormals->SplittingOff();
  marginalNormals->SetInput(surface);
  marginalNormals->Update();
  surface = marginalNormals->GetOutput();
  surface->Modified();

  // Actually do the writing.
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(marginalMeanOutputFile);
  writer->Write();

  // Now repeat for the 'dumb' mean.
  for (i = 0; i < surface->GetNumberOfPoints(); i++)
  {
    coord[0] = MeanShape_2(3 * i    );
    coord[1] = MeanShape_2(3 * i + 1);
    coord[2] = MeanShape_2(3 * i + 2);
    surface->GetPoints()->SetPoint(i, coord);
  }
  surface->Modified();

  vtkPolyDataNormals *dumbNormals = vtkPolyDataNormals::New();
  dumbNormals->SplittingOff();
  dumbNormals->SetInput(surface);
  dumbNormals->Update();
  surface = dumbNormals->GetOutput();
  surface->Modified();

  vtkPolyDataWriter *writer2 = vtkPolyDataWriter::New();
  writer2->SetInput(surface);
  writer2->SetFileName(dumbMeanOutputFile);
  writer2->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library"
       << endl;
}
#endif
