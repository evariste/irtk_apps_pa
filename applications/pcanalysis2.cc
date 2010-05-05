
#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <stdio.h>
#include <string.h>
#include <irtkImage.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>

// Minimum norm of an eigenvector
#define MIN_NORM 0.01

// vtk includes
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

vtkPointSet *data;
vtkPoints  *points;
vtkDataArray *vectors;

char* componentsFile = NULL;
char* evecsFile = NULL;
char* evalsFile = NULL;
char* meanFile  = NULL;

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
	writer->SetFileTypeToBinary();
	writer->Update();
      } else {
	cerr << "Unknown VTK data type" << endl;
	exit(1);
      }
    }
  }
}

int main(int argc, char **argv)
{
  int i, j, k, m, iNoOfLandmarks, iNoOfDatasets;
  irtkMatrix T, Eigenvector;
  irtkVector Da, Db, MeanShape, Eigenvalues;
  double p[3], v[3];
  int maxModes = 10;

  if ((argc < 5)){
    cout << "usage: pcanalysis landmarks1 ... landmarks_N landmarks_mean "
	 << "eigen_values eigen_vectors maxModes componentMatrix \n";
    return 1;
  }

  componentsFile = argv[argc-1];
  maxModes = atoi(argv[argc-2]);
  evecsFile = argv[argc-3];
  evalsFile = argv[argc-4];
  meanFile  = argv[argc-5];

  cout << "Using " << maxModes << " modes\n";

  iNoOfDatasets  = argc-6;
  iNoOfLandmarks = 0;
  cerr << " There are " << iNoOfDatasets << " data sets for training"
       << endl;

  // Read all landmarks in matrix M
  for (i = 0; i < iNoOfDatasets; i++){
    cerr << " Including landmarks in " << argv[i+1] << endl;
    data    = read(argv[i+1]);
    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks = points->GetNumberOfPoints();
//       M.Initialize(3*iNoOfLandmarks,iNoOfDatasets);
      Da.Initialize(3*iNoOfLandmarks);
      Db.Initialize(3*iNoOfLandmarks);

      MeanShape.Initialize(3*iNoOfLandmarks);
      for (j = 0; j < 3*iNoOfLandmarks; j++){
        MeanShape(j) = 0.0;
      }
      cerr << " There are " << iNoOfDatasets 
           << " datasets with "
           << points->GetNumberOfPoints() 
           << " landmarks." << endl;
    } else {
      if ((points->GetNumberOfPoints()) != iNoOfLandmarks){
        cerr << "Datasets must contain the same number of landmarks" << endl;
        exit(1);          
      }
    }
    
    for (j = 0; j < iNoOfLandmarks; j++){
      if (vectors == NULL){
	points->GetPoint(j,p);
	MeanShape(3*j)   += p[0]; 
	MeanShape(3*j+1) += p[1]; 
	MeanShape(3*j+2) += p[2]; 
      } else {
	vectors->GetTuple(j,v);
	MeanShape(3*j)   += v[0]; 
	MeanShape(3*j+1) += v[1]; 
	MeanShape(3*j+2) += v[2]; 
      }
    }
  }

  MeanShape = MeanShape*(1.0/iNoOfDatasets);
  cout << "Calculated mean." << endl;
  
  // Form matrix T = (1/iNoOfDatasets) M^T * M which has
  // dimensions iNoOfDatasets x iNoOfDatasets after subtracting
  // the mean from each dataset.

  T.Initialize(iNoOfDatasets,iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; ++i){

    data    = read(argv[i+1]);
    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    for (k = 0; k < iNoOfLandmarks; ++k){
      if (vectors == NULL){
	points->GetPoint(k,p);
 	Da(3*k)   = p[0] - MeanShape(3*k); 
	Da(3*k+1) = p[1] - MeanShape(3*k+1); 
	Da(3*k+2) = p[2] - MeanShape(3*k+2); 
      } else {
	vectors->GetTuple(k,v);
 	Da(3*k)   = v[0] - MeanShape(3*k); 
	Da(3*k+1) = v[1] - MeanShape(3*k+1); 
	Da(3*k+2) = v[2] - MeanShape(3*k+2); 
      }
    }

    data->Delete();

    for (j = i; j < iNoOfDatasets; ++j){

      data    = read(argv[j + 1]);
      points  = data->GetPoints();
      vectors = data->GetPointData()->GetVectors();

      for (k = 0; k < iNoOfLandmarks; ++k){
        if (vectors == NULL){
          points->GetPoint(k,p);
          Db(3*k)   = p[0] - MeanShape(3*k);
          Db(3*k+1) = p[1] - MeanShape(3*k+1);
          Db(3*k+2) = p[2] - MeanShape(3*k+2);
        } else {
          vectors->GetTuple(k,v);
          Db(3*k)   = v[0] - MeanShape(3*k);
          Db(3*k+1) = v[1] - MeanShape(3*k+1);
          Db(3*k+2) = v[2] - MeanShape(3*k+2);
        }
      }

      data->Delete();

      // Covariance matrix entry
      T(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks; k++){
        T(i, j) += Da(k) * Db(k);
      }

      // Symmetric covariance matrix
      T(j, i) = T(i, j);

    }
    
  }

  T /= iNoOfDatasets;
    
  // Compute the eigen decomposition
  irtkEigenAnalysis ea(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++){
    for(int j=0; j < iNoOfDatasets; j++){
      ea.Matrix(i,j) = T(i,j);  
    }
  }
  ea.DecrSortEigenStuff();

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M*M^T.
  Eigenvector.Initialize(3*iNoOfLandmarks, maxModes);
  Eigenvector *= 0.0;

  for (k = 0; k < iNoOfDatasets; k++){
    data    = read(argv[k+1]);
    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    for (m = 0; m < iNoOfLandmarks; ++m){
      if (vectors == NULL){
        points->GetPoint(m,p);
        Da(3*m)   = p[0] - MeanShape(3*m);
        Da(3*m+1) = p[1] - MeanShape(3*m+1);
        Da(3*m+2) = p[2] - MeanShape(3*m+2);
      } else {
        vectors->GetTuple(m,v);
        Da(3*m)   = v[0] - MeanShape(3*m);
        Da(3*m+1) = v[1] - MeanShape(3*m+1);
        Da(3*m+2) = v[2] - MeanShape(3*m+2);
      }
    }

    data->Delete();

    for (i = 0; i < 3*iNoOfLandmarks; ++i){
      for (j = 0; j < maxModes; ++j){
        Eigenvector(i, j) += Da(i) * ea.Eigenvector(k, j);
      }
    }
  }

  Eigenvalues.Initialize(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets ; i++){
    Eigenvalues(i) = ea.Eigenvalue(i);
  }
  
  float fTotalVar   = 0;
  float fCummulated = 0;
  for (i = 0; i < iNoOfDatasets; i++){
    fTotalVar += ea.Eigenvalue(i);
  }
  for (i = 0; i < maxModes; i++){
    cout<< "Mode:" << i << " Eigenvalue:"<<ea.Eigenvalue(i)<<endl;

    fCummulated += 100 * ea.Eigenvalue(i) / fTotalVar;
    cout << " Mode [" << i << "] explains " 
         << 100 * ea.Eigenvalue(i) / fTotalVar
         << " % (" 
         << fCummulated 
         << " %) of shape variance." << endl;
  }

  // Normalize eigenvectors
  for (i = 0; i < maxModes ; i++){
    float fNorm = 0.0;
    for (j = 0; j < 3*iNoOfLandmarks; j++){
      fNorm += Eigenvector(j,i) * Eigenvector(j,i);
    }
    fNorm = sqrt(fNorm);
    if (100 * ea.Eigenvalue(i) / fTotalVar > MIN_NORM){
      for (j = 0; j < 3*iNoOfLandmarks; j++){
	Eigenvector(j,i) /= fNorm;
      }
    } else {
      for (j = 0; j < 3*iNoOfLandmarks; j++){
	Eigenvector(j,i) = 0;
      }
      Eigenvalues(i) = 0;
    }
  }

  cout << "Eigenvectors normalised" << endl;

  // Find components matrix.
  irtkMatrix comps;
  comps.Initialize(maxModes, iNoOfDatasets);

  for (k = 0; k < iNoOfDatasets; k++){
    data    = read(argv[k+1]);
    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    for (m = 0; m < iNoOfLandmarks; ++m){
      if (vectors == NULL){
        points->GetPoint(m,p);
        Da(3*m)   = p[0] - MeanShape(3*m);
        Da(3*m+1) = p[1] - MeanShape(3*m+1);
        Da(3*m+2) = p[2] - MeanShape(3*m+2);
      } else {
        vectors->GetTuple(m,v);
        Da(3*m)   = v[0] - MeanShape(3*m);
        Da(3*m+1) = v[1] - MeanShape(3*m+1);
        Da(3*m+2) = v[2] - MeanShape(3*m+2);
      }
    }

    data->Delete();

    for (j = 0; j < maxModes; ++j){
      comps(j, k) = 0;
      for (i = 0; i < 3*iNoOfLandmarks; ++i){
        comps(j, k) += Eigenvector(i,j) * Da(i);
      }
    }

  }


  cerr << " Writing eigenvalues to " << evalsFile << endl;
  Eigenvalues.Write(evalsFile);
    
  cerr << " Writing eigenshapes to " << evecsFile << endl;
  Eigenvector.Write(evecsFile);

  cout << "Writing components to " << componentsFile << endl;
  comps.Write(componentsFile);
  
  
  // Read the first dataset for its structure.
  data    = read(argv[1]);
  points  = data->GetPoints();
  vectors = data->GetPointData()->GetVectors();

  // Convert mean back to VTK
  for (j = 0; j < iNoOfLandmarks; j++){
    double p[3], v[3];
    if (vectors == NULL){
      p[0] = MeanShape(3*j);
      p[1] = MeanShape(3*j+1); 
      p[2] = MeanShape(3*j+2); 
      points->SetPoint(j,p);
    } else {
      v[0] = MeanShape(3*j);
      v[1] = MeanShape(3*j+1); 
      v[2] = MeanShape(3*j+2); 
      vectors->SetTuple(j,v);
    }	 
  }
  
  data->SetPoints(points);
  cerr << " Writing mean shape to " << meanFile << endl;  
  write(meanFile, data);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library" 
       << endl;
}
#endif
