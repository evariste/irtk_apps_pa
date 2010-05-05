#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <stdio.h>
#include <string.h>
#include <irtkImage.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>

// Minimum norm of an eigenvector
#define MIN_NORM 0.01

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkStructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

vtkPointSet *data;
vtkPoints  *points;
vtkDataArray *vectors;

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

void usage(){
  cout << "usage: pcanalysis-variances-only N landmarks_1 ... landmarks_N <options>" << "\n";
  cout << "  options:" << endl;
  cout << "  -q               quiet output, don't list details for all modes."  << endl;
  cout << "  -percent value   Return the number of modes that explain " << endl;
  cout << "                   <value>% of the variance" << endl;
  cout << "  -csv             csv output, forces -q and reports 'datasets,landmarks,totalVar,meanVar,geoMean'" << endl;
  cout << "                   this is appended with 'pcaError,percent,modesRequired' if '-percent' is specified." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, iNoOfLandmarks, iNoOfDatasets;
  int verbose = True;
  int ok;
  int csvMode = False;
  double norm;
  int modesRequired;

  int percent = -1;
  double retrievedComp, squaredDiff;
  double totalVar, totalLogVar, accumulated, totalLogCount;

  irtkMatrix M, T, Eigenvector;
  irtkVector MeanShape, Eigenvalues;
  irtkMatrix data_reduced;

  if ((argc < 4)){
    usage();
  }

  iNoOfDatasets = atoi(argv[1]);
  argv++;
  argc--;

  iNoOfLandmarks = 0;

  modesRequired = iNoOfDatasets;

  // Read all landmarks in matrix M
  for (i = 0; i < iNoOfDatasets; ++i){
    data    = read(argv[1]);
    argv++;
    argc--;

    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks = points->GetNumberOfPoints();
      M.Initialize(3*iNoOfLandmarks,iNoOfDatasets);
      MeanShape.Initialize(3*iNoOfLandmarks);
      for (j = 0; j < 3*iNoOfLandmarks; ++j){
        MeanShape(j) = 0.0;
      }

    } else {
      if ((points->GetNumberOfPoints()) != iNoOfLandmarks){
        cerr << "Datasets must contain the same number of landmarks" << endl;
        exit(1);          
      }
    }
    
    for (j = 0; j < iNoOfLandmarks; ++j){
      double p[3], v[3];
      if (vectors == NULL){
	points->GetPoint(j, p);
	M(3*j,i)   = p[0];
	M(3*j+1,i) = p[1];
	M(3*j+2,i) = p[2];
	MeanShape(3*j)   += p[0]; 
	MeanShape(3*j+1) += p[1]; 
	MeanShape(3*j+2) += p[2]; 
      } else {
	vectors->GetTuple(j, v);
	M(3*j,i)   = v[0];
	M(3*j+1,i) = v[1];
	M(3*j+2,i) = v[2];
	MeanShape(3*j)   += v[0]; 
	MeanShape(3*j+1) += v[1]; 
	MeanShape(3*j+2) += v[2]; 
      }
    }
  }

  MeanShape = MeanShape*(1.0/iNoOfDatasets);
  
  // Read optional arguments.
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      verbose = False;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-percent") == 0)){
      argc--;
      argv++;
      percent = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-csv") == 0)){
      argc--;
      argv++;
      csvMode = True;
      ok = True;
    }

    if (ok == False){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (csvMode == True){
    verbose = False;
  }

  if (percent > -1){
    if (percent > 100){
      cerr << "Percentage must be between 0 and 100" << endl;
      exit (1);
    }
  }

  if (verbose == True){
    cerr << " There are " << iNoOfDatasets;
    cout << " datasets with " << iNoOfLandmarks << " landmarks." << endl;
  }

  // Subtract the mean  
  for (i = 0; i < iNoOfDatasets; ++i){
    for (j = 0; j < iNoOfLandmarks; ++j){
      M(3 * j    , i) -= MeanShape(3*j);
      M(3 * j + 1, i) -= MeanShape(3*j+1);
      M(3 * j + 2, i) -= MeanShape(3*j+2);
    }
  }

  // Form matrix T = (1/iNoOfDatasets) M^T * M which has dimensions
  // iNoOfDatasets x iNoOfDatasets
  T.Initialize(iNoOfDatasets,iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      T(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks; ++k)
        T(i, j) += M(k, i) * M(k, j);
    }
  }

  T /= iNoOfDatasets;
    
  // Compute the eigen decomposition
  irtkEigenAnalysis ea(iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets; ++i){
    for(j = 0; j < iNoOfDatasets; ++j){
      ea.Matrix(i, j) = T(i, j);  
    }
  }
  ea.DecrSortEigenStuff();

  ////////////////////////////////////////////////////////

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M*M^T.
  Eigenvector.Initialize(3 * iNoOfLandmarks, iNoOfDatasets);
  Eigenvector *= 0.0;

  for (i = 0; i < 3 * iNoOfLandmarks ; ++i){
    for (j = 0; j < iNoOfDatasets; ++j){
      for (k = 0; k < iNoOfDatasets; ++k){
        Eigenvector(i, j) += M(i, k) * ea.Eigenvector(k, j);
      }
    }
  }

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets ; ++i){
    norm = 0.0;

    for (j = 0; j < 3 * iNoOfLandmarks; ++j){
      norm += Eigenvector(j, i) * Eigenvector(j, i);
    }
    norm = sqrt(norm);

    for (j = 0; j < 3 * iNoOfLandmarks; ++j){
      Eigenvector(j, i) /= norm;
    }
  }

  /////////////////////////////////////////////////////

  Eigenvalues.Initialize(iNoOfDatasets);

  for (i = 0; i < iNoOfDatasets ; ++i){
    Eigenvalues(i) = ea.Eigenvalue(i);
  }
 
  totalVar = totalLogVar = accumulated = totalLogCount = 0.0;

  for (i = 0; i < iNoOfDatasets; ++i){

    totalVar += ea.Eigenvalue(i);

    if (ea.Eigenvalue(i) > 0){
      totalLogVar += log(ea.Eigenvalue(i));
      totalLogCount++;
    }

  }

  if (verbose == True){
    cout << "Mode\tEigenvalue\texplains\tcumulative " << endl;
  }

  for (i = 0; i < iNoOfDatasets; ++i){
    if (verbose == True){
      cout << i << "\t" <<ea.Eigenvalue(i);
      accumulated += 100 * ea.Eigenvalue(i) / totalVar;
      cout << "\t" << 100 * ea.Eigenvalue(i) / totalVar;
      cout << "\t" << accumulated << endl;
    }
  }

  if (csvMode == True){
    cout << iNoOfDatasets << "," << iNoOfLandmarks << "," << totalVar << ",";
    cout << totalVar / iNoOfDatasets << ",";
    cout << exp (totalLogVar / totalLogCount);
  } else {
    cout << "Total variance = " << totalVar << endl;
    cout << "Mean  variance = " << totalVar / iNoOfDatasets << endl;
    cout << "Geomean        = " << exp (totalLogVar / totalLogCount) << endl;
  }

  if (percent > -1){

    accumulated = 0.0;
    i = 0;

    while (accumulated < percent && i < iNoOfDatasets){
      accumulated += 100 * ea.Eigenvalue(i) / totalVar;
      ++i;
    }

    modesRequired = i;

    // Now convert x-u data into pca coordinates using only first R eigenvectors.
    data_reduced.Initialize(modesRequired, iNoOfDatasets);

    for (i = 0; i <  modesRequired; ++i){
      for (j = 0; j < iNoOfDatasets; ++j){
        data_reduced(i, j) = 0.0;
        for (k = 0; k < 3 * iNoOfLandmarks; ++k)
          data_reduced(i, j) += Eigenvector(k, i) * M(k, j);
      }
    }

    // data_reduced(i, j) is the component of the jth column of M (jth
    // shape) w.r.t the ith eigenvector.
    // I.e. each column of data_reduced are the components of one shape.

    // Following are checks to see whether projections into eigenspaces worked.
    // Check data transformation is working by transforming reduced data back
    // and comparing with original data.

    squaredDiff = 0.0;

    for (i = 0; i < 3 * iNoOfLandmarks; ++i){
      for (j = 0; j < iNoOfDatasets; ++j){

        retrievedComp = 0.0;

        for (k = 0; k < modesRequired; ++k){
          retrievedComp += Eigenvector(i, k) * data_reduced(k, j);
        }

        squaredDiff += ( retrievedComp - M(i,j) ) * ( retrievedComp - M(i, j) );

      }
    }

    if (csvMode == True){
      cout << "," << squaredDiff / (iNoOfLandmarks * iNoOfDatasets);
      cout << "," << percent << "," << modesRequired;
    } else {
      cout << "PCA reconstruction error = " << squaredDiff / (iNoOfLandmarks * iNoOfDatasets);
      cout << modesRequired  << " modes required to explain at least " << percent << "% of the variance" << endl;
      cout << "Mean squared distance using truncated modes =  " << squaredDiff / (iNoOfLandmarks * iNoOfDatasets) << endl;
    }

  } // if (percent > -1)

  cout << endl;

  ////////////////////////////////////////////////////
  // Debugging

  // For some reason, some matrices do not converge when the eigenvectors
  // are being calculated. The above code uses the QL algorithm.  Also
  // checked using SVD below, because the covariance matrices are symmetric
  // and positive definite (hopefully), the SVD should give the eigen
  // stuff.  Matrices that do not converge using the NR QL method above
  // also fail to converge for the SVD method below.  The eigen values and
  // vectors still appear reasonable however.

  // The above seemed to be fixed by increasing the number of iterations in
  // the irtkEigenanalysis::QL routine.

//    Eigenvector.Write("evm.mat");
//    data_reduced.Write("truncpcacomps.mat");
//    T.Write("t.mat");
//    M.Write("m.mat");

//   irtkMatrix u, v, d;
//   irtkVector w;

//   u.Initialize(iNoOfDatasets, iNoOfDatasets);
//   v.Initialize(iNoOfDatasets, iNoOfDatasets);
//   w.Initialize(iNoOfDatasets);

//   cout << "Initialised matrices." << endl;

//   T.SVD(v, w, u);
//   // m = v * d * (~u);
//   // where d is a diagonal matrix with entries from w.

//   cout << "svd done." << endl;

//   for (i = 0; i < iNoOfDatasets; ++i){
//     cout << ea.Eigenvalue(i) << "\t" << w(i) << endl;
//   }

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library" 
       << endl;
}
#endif


////////////////////////////////////////////////////////
// EXAMPLE OF MATLAB CODE THAT DOES SOMETHING SIMILAR TO THE ABOVE.  MAINLY
// TO CHECK THE SQUARED DIFFERENCES WERE BEING CALCULATED CORRECTLY.
////////////////////////////////////////////////////////

/*

    noOfDatasets = 178;
    noOfLandmarks = 1630;

    %Read the data from the c++ exe.
    m = readItkMatrix('m.mat');
    %Read the covariance matrix from the c++ exe.
    t = readItkMatrix('t.mat');

    % Do a PCA by finding the eigenvalues and vectors of the cov matrix.
    [ut, st, vt] = svd(t);

    % Convert the eigenvectors of the covariance matrix (small dimensional)
    % into eigenvectors of the original data (high dimensional).
    evm = m * vt;

    % Normalise the eigenvectors.
    for i = 1:noOfDatasets
    evm(:,i) = evm(:, i) / norm(evm(:, i));
    end

    % What are the PCA components of the original data w.r.t the basis
    % defined by the eigenvectors?
    pcacomps = m.' * evm;

    % Let's truncate some modes and try and reconstruct the data using only
    % the modes that are kept.
    modeskept = 5;

    retrieved = evm(:, 1:modeskept) * (pcacomps(:, 1:modeskept)).';

    % Compare with the original data.
    diff = m - retrieved;
    diffsq = diff .* diff;

    % Average squared difference for a landmark across all shapes.
    sum(sum(diffsq)) / noOfDatasets / noOfLandmarks

    %%%%%%%%%%%%%%%%%%%%%%%
    % Miscellaneous checks.
    %%%%%%%%%%%%%%%%%%%%%%%
    disp ('Checking that the covariance matrix produced by the c++ exe is accurate.');
    diff = (m.' * m ) / noOfDatasets - t;
    max(max(abs(diff)))

    % Because the covariance matrix is symmetric and positive definite, the
    % right and left singular vectors should be equal to each other and equal
    % to the eigenvectors of T
    diff = ut - vt;
    disp ('Mean diff between right and left singular vectors.')
    mean(mean(abs(diff)))
    %imagesc(diff)

    % Checking that the columns of vt represent the eigenvectors of T.
    et = diag(st);
    diff = t * vt(:, 1) - et(1) * vt(:, 1);
    disp('checking that the singular vectors of T are eigenvectors.');
    max(abs(diff))

*/



