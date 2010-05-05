#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <fstream>
#include <stdio.h>
#include <string.h>
#include <irtkImage.h>
// #include <irtkEigenAnalysis.h>
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
char *templateFileName = NULL;
char **outputFilenames = NULL;
char **inputFilenames  = NULL;
char *meanFilename     = NULL;

int alignVersion = 1;

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
	writer->SetFileTypeToASCII();
	writer->Write();
      } else {
	cerr << "Unknown VTK data type" << endl;
	exit(1);
      }
    }
  }
}

void print_matrix(irtkMatrix M)
{
 int i,j;
 cerr << "Matrix start" << endl;
 for (i=0;i<M.Rows();i++){
 	for (j=0;j<M.Cols();j++){
		if (j!=M.Cols()-1) cerr << M(i,j) << " ";
		else cerr << M(i,j) <<endl;
	}
 }  
 cerr << "Matrix end" << endl;
}


// Returns result of aligning M_1 to M_2.
// Reference argument transformation will contain the homogeneous transformation matrix.
irtkMatrix align2(irtkMatrix & M_2, irtkMatrix & M_1, int iNoOfLandmarks, irtkMatrix & transformationMatrix)
{
  // calculate centroids
  irtkMatrix centroid_1, centroid_2;
  centroid_1.Initialize(3, 1);
  centroid_2.Initialize(3, 1);

  int j;

  for (j = 0; j < iNoOfLandmarks; j++){
    centroid_1(0, 0) += M_1(3 * j, 0);
    centroid_1(1, 0) += M_1(3 * j + 1, 0);
    centroid_1(2, 0) += M_1(3 * j + 2, 0); 
  }

  for (j = 0; j < iNoOfLandmarks; j++){
    centroid_2(0, 0) += M_2(3 * j, 0);
    centroid_2(1, 0) += M_2(3 * j + 1, 0);
    centroid_2(2, 0) += M_2(3 * j + 2, 0); 
  }

  centroid_1 *= (1.0/iNoOfLandmarks);
  centroid_2 *= (1.0/iNoOfLandmarks);

  irtkVector M_1_centred, M_2_centred;

  M_1_centred.Initialize(3 * iNoOfLandmarks);
  M_2_centred.Initialize(3 * iNoOfLandmarks);

  for (j = 0; j < iNoOfLandmarks; j++){
    M_1_centred(3 * j)     = M_1(3 * j    , 0) - centroid_1(0, 0);
    M_1_centred(3 * j + 1) = M_1(3 * j + 1, 0) - centroid_1(1, 0);
    M_1_centred(3 * j + 2) = M_1(3 * j + 2, 0) - centroid_1(2, 0);

    M_2_centred(3 * j)     = M_2(3 * j    , 0) - centroid_2(0, 0);
    M_2_centred(3 * j + 1) = M_2(3 * j + 1, 0) - centroid_2(1, 0);
    M_2_centred(3 * j + 2) = M_2(3 * j + 2, 0) - centroid_2(2, 0);
  }

  double norm_1 = 0.0;
  double norm_2 = 0.0;

  for (j = 0; j < 3 * iNoOfLandmarks; ++j){
    norm_1 += M_1_centred(j) * M_1_centred(j);
    norm_2 += M_2_centred(j) * M_2_centred(j);
  }

  norm_1 = sqrt(norm_1);
  norm_2 = sqrt(norm_2);

  // Normalise the scales of the datasets.
  for (j = 0; j < 3 * iNoOfLandmarks; ++j){
    M_1_centred(j) = M_1_centred(j) / norm_1;
    M_2_centred(j) = M_2_centred(j) / norm_2;
  }

  // Treat M_1_centred and M_2_centred as Z_1 and Z_2 as per Mardia & Dryden.
  // calculate Z_2^T * Z_1, call it m.

  irtkMatrix m;
  m.Initialize(3, 3);

  for (j = 0; j < iNoOfLandmarks; j++){
    m(0, 0) += M_2_centred(3 * j    ) * M_1_centred(3 * j    );
    m(0, 1) += M_2_centred(3 * j    ) * M_1_centred(3 * j + 1);
    m(0, 2) += M_2_centred(3 * j    ) * M_1_centred(3 * j + 2);

    m(1, 0) += M_2_centred(3 * j + 1) * M_1_centred(3 * j    );
    m(1, 1) += M_2_centred(3 * j + 1) * M_1_centred(3 * j + 1);
    m(1, 2) += M_2_centred(3 * j + 1) * M_1_centred(3 * j + 2);

    m(2, 0) += M_2_centred(3 * j + 2) * M_1_centred(3 * j    );
    m(2, 1) += M_2_centred(3 * j + 2) * M_1_centred(3 * j + 1);
    m(2, 2) += M_2_centred(3 * j + 2) * M_1_centred(3 * j + 2);
  }


  irtkMatrix u, v, d;
  irtkVector w;

  d.Initialize(3, 3);
  u.Initialize(3, 3);
  v.Initialize(3, 3);
  w.Initialize(3);

  m.SVD(v, w, u);
  // m = v * d * (~u);

  d(0, 0) = w(0);
  d(1, 1) = w(1);
  d(2, 2) = w(2);

  double scale, beta;

  beta = w(0) + w(1) + w(2);
  scale = norm_2 * beta / norm_1;
  // Ignoring beta seems to give better results.
  scale = norm_2 / norm_1;

  irtkMatrix rot;
  rot.Initialize(3, 3);
  rot = u * (~v);

  // now calculate translation t
  irtkMatrix t;
  t.Initialize(3, 1);
  t = centroid_2 - (~rot) * centroid_1 * scale;

  // now transform shape 2 into coordinate system of shape 1
  irtkMatrix transformed_M_1;
  transformed_M_1.Initialize(3 * iNoOfLandmarks, 1);

  for (j = 0; j < iNoOfLandmarks; j++){
    transformed_M_1(3 * j, 0)     = M_1(3 * j, 0) * rot(0, 0) + M_1(3 * j + 1, 0) * rot(1, 0) + M_1(3 * j + 2, 0) * rot(2, 0);
    transformed_M_1(3 * j + 1, 0) = M_1(3 * j, 0) * rot(0, 1) + M_1(3 * j + 1, 0) * rot(1, 1) + M_1(3 * j + 2, 0) * rot(2, 1);
    transformed_M_1(3 * j + 2, 0) = M_1(3 * j, 0) * rot(0, 2) + M_1(3 * j + 1, 0) * rot(1, 2) + M_1(3 * j + 2, 0) * rot(2, 2);
  }

  transformed_M_1 *= scale;

  for (j = 0; j < iNoOfLandmarks; j++){
    transformed_M_1(3 * j, 0)     += t(0, 0);
    transformed_M_1(3 * j + 1, 0) += t(1, 0);
    transformed_M_1(3 * j + 2, 0) += t(2, 0);
  }


  // Store the transformation matrix if it is needed.
  transformationMatrix *= 0;
  // Put the rotation matrix into the upper left corner.
  rot.Transpose();
  transformationMatrix(rot, 0, 0);
  transformationMatrix *= scale;
  // Put the translation in the right column
  transformationMatrix(t, 0, 3);
  // Set the bottom right corner.
  transformationMatrix(3, 3) = 1.0;

  return transformed_M_1;

}

// Returns result of aligning M_2 to M_1.
// Reference argument transformation will contain the homogeneous transformation matrix.
irtkMatrix align1(irtkMatrix & M_1, irtkMatrix & M_2, int iNoOfLandmarks, irtkMatrix & transformationMatrix)
{
  // calculate centroids
  irtkMatrix centroid_1, centroid_2;
  centroid_1.Initialize(3, 1);
  centroid_2.Initialize(3, 1);

  int j;

  for (j = 0; j < iNoOfLandmarks; j++){
    centroid_1(0, 0) += M_1(3 * j, 0);
    centroid_1(1, 0) += M_1(3 * j + 1, 0);
    centroid_1(2, 0) += M_1(3 * j + 2, 0); 
  }

  for (j = 0; j < iNoOfLandmarks; j++){
    centroid_2(0, 0) += M_2(3 * j, 0);
    centroid_2(1, 0) += M_2(3 * j + 1, 0);
    centroid_2(2, 0) += M_2(3 * j + 2, 0); 
  }

  centroid_1 *= (1.0/iNoOfLandmarks);
  centroid_2 *= (1.0/iNoOfLandmarks);

  irtkVector M_1_centred, M_2_centred;

  M_1_centred.Initialize(3 * iNoOfLandmarks);
  M_2_centred.Initialize(3 * iNoOfLandmarks);

  for (j = 0; j < iNoOfLandmarks; j++){
    M_1_centred(3 * j)     = M_1(3 * j    , 0) - centroid_1(0, 0);
    M_1_centred(3 * j + 1) = M_1(3 * j + 1, 0) - centroid_1(1, 0);
    M_1_centred(3 * j + 2) = M_1(3 * j + 2, 0) - centroid_1(2, 0);

    M_2_centred(3 * j)     = M_2(3 * j    , 0) - centroid_2(0, 0);
    M_2_centred(3 * j + 1) = M_2(3 * j + 1, 0) - centroid_2(1, 0);
    M_2_centred(3 * j + 2) = M_2(3 * j + 2, 0) - centroid_2(2, 0);
  }

  // have centroid centred data

  // now calculate sums for matrix N, shape 1 is right coord system, shape
  // 2 is left coord system, doing left to right

  float sxx, sxy, sxz, syx, syy, syz, szx, szy, szz;
  sxx = sxy = sxz = syx = syy = syz = szx = szy = szz = 0;

  for (j = 0; j < iNoOfLandmarks; j++){
    sxx += M_2_centred(3 * j)     * M_1_centred(3 * j);
    sxy += M_2_centred(3 * j)     * M_1_centred(3 * j + 1);
    sxz += M_2_centred(3 * j)     * M_1_centred(3 * j + 2);

    syx += M_2_centred(3 * j + 1) * M_1_centred(3 * j);
    syy += M_2_centred(3 * j + 1) * M_1_centred(3 * j + 1);
    syz += M_2_centred(3 * j + 1) * M_1_centred(3 * j + 2);

    szx += M_2_centred(3 * j + 2) * M_1_centred(3 * j);
    szy += M_2_centred(3 * j + 2) * M_1_centred(3 * j + 1);
    szz += M_2_centred(3 * j + 2) * M_1_centred(3 * j + 2);
  }

  float **N;

  N = matrix(1, 4, 1, 4);
  N[1][1] =  sxx + syy + szz;
  N[2][2] =  sxx - syy - szz;
  N[3][3] = -sxx + syy - szz;
  N[4][4] = -sxx - syy + szz;

  N[1][2] = N[2][1] = syz - szy;
  N[1][3] = N[3][1] = szx - sxz;
  N[2][3] = N[3][2] = sxy + syx;
  N[1][4] = N[4][1] = sxy - syx;
  N[2][4] = N[4][2] = szx + sxz;
  N[3][4] = N[4][3] = syz + szy;

  int nrot;
  float *er;
  float **ev;

  er = vector(1, 4);
  ev = matrix(1, 4, 1, 4);

  jacobi(N, 4, er, ev,&nrot);
  // sort by descending eigenvalue
  eigsrt(er, ev, 4); 

  //1st eigenvector is quaternion of the rotation
  irtkVector quat;
  quat.Initialize(4);
  quat(0) = ev[1][1];
  quat(1) = ev[2][1];
  quat(2) = ev[3][1];
  quat(3) = ev[4][1];
  // calculate rotation matrix

  irtkMatrix Rot;
  Rot.Initialize(3, 3);
  Rot(0, 0) = quat(0) * quat(0) + quat(1) * quat(1) - quat(2) * quat(2) - quat(3) * quat(3);
  Rot(1, 1) = quat(0) * quat(0) - quat(1) * quat(1) + quat(2) * quat(2) - quat(3) * quat(3);
  Rot(2, 2) = quat(0) * quat(0) - quat(1) * quat(1) - quat(2) * quat(2) + quat(3) * quat(3);

  Rot(0, 1) = 2 * (quat(1) * quat(2) - quat(0) * quat(3));
  Rot(1, 0) = 2 * (quat(1) * quat(2) + quat(0) * quat(3));
  Rot(0, 2) = 2 * (quat(1) * quat(3) + quat(0) * quat(2));
  Rot(2, 0) = 2 * (quat(1) * quat(3) - quat(0) * quat(2));
  Rot(1, 2) = 2 * (quat(2) * quat(3) - quat(0) * quat(1));
  Rot(2, 1) = 2 * (quat(2) * quat(3) + quat(0) * quat(1));

  // calculate scale, symmetric scaling is used

  float scale, shape1dev, shape2dev;
  shape1dev = shape2dev = 0.0;
  for (j = 0; j < iNoOfLandmarks; j++){
        shape1dev += M_1_centred(3 * j)     * M_1_centred(3 * j);
        shape1dev += M_1_centred(3 * j + 1) * M_1_centred(3 * j + 1);
        shape1dev += M_1_centred(3 * j + 2) * M_1_centred(3 * j + 2);
        shape2dev += M_2_centred(3 * j)     * M_2_centred(3 * j);
        shape2dev += M_2_centred(3 * j + 1) * M_2_centred(3 * j + 1);
        shape2dev += M_2_centred(3 * j + 2) * M_2_centred(3 * j + 2);
  }

  scale = sqrt(shape1dev/shape2dev);

  // now calculate translation t
  irtkMatrix t;
  t.Initialize(3, 1);
  t = centroid_1 - Rot * centroid_2 * scale;

  // now transform shape 2 into coordinate system of shape 1

  irtkMatrix transformed_M_2;
  transformed_M_2.Initialize(3 * iNoOfLandmarks, 1);

  for (j = 0; j < iNoOfLandmarks; j++){
   transformed_M_2(3 * j, 0)     = Rot(0, 0) * M_2(3 * j, 0) + Rot(0, 1) * M_2(3 * j + 1, 0) + Rot(0, 2) * M_2(3 * j + 2, 0);
   transformed_M_2(3 * j + 1, 0) = Rot(1, 0) * M_2(3 * j, 0) + Rot(1, 1) * M_2(3 * j + 1, 0) + Rot(1, 2) * M_2(3 * j + 2, 0);
   transformed_M_2(3 * j + 2, 0) = Rot(2, 0) * M_2(3 * j, 0) + Rot(2, 1) * M_2(3 * j + 1, 0) + Rot(2, 2) * M_2(3 * j + 2, 0);
  }

  transformed_M_2 *= scale;

  for (j = 0; j < iNoOfLandmarks; j++){
   transformed_M_2(3 * j, 0)     = transformed_M_2(3 * j, 0)     + t(0, 0);
   transformed_M_2(3 * j + 1, 0) = transformed_M_2(3 * j + 1, 0) + t(1, 0);
   transformed_M_2(3 * j + 2, 0) = transformed_M_2(3 * j + 2, 0) + t(2, 0);
  }

  // Store the transformation matrix if it is needed.
  transformationMatrix *= 0;
  // Put the rotation matrix into the upper left corner.
  transformationMatrix(Rot, 0, 0);
  transformationMatrix *= scale;
  // Put the translation in the right column
  transformationMatrix(t, 0, 3);
  // Set the bottom right corner.
  transformationMatrix(3, 3) = 1.0;

  return transformed_M_2;
}

irtkMatrix align(irtkMatrix & M_1, irtkMatrix & M_2, int iNoOfLandmarks, irtkMatrix & transformationMatrix)
{

  if (alignVersion == 1){
    return align1(M_1, M_2, iNoOfLandmarks, transformationMatrix);
  } else if (alignVersion == 2) {
    return align2(M_1, M_2, iNoOfLandmarks, transformationMatrix);
  } else {
    cerr << "non-valid align version : " << alignVersion << endl;
    exit(1);
  }

}


void usage(){
  cerr << "Usage: procrustes_align N [inSurf_1] ... [inSurf_N] [outSurf_1] ... [outSurf_N] <options>" << endl;
  cerr << "Procrustes alignment of inSurf_1 ... inSurf_N.  Write results to outSurf_1 ... outSurf_N" << endl;
  cerr << "N is the number of datasets." << endl;
  cerr << "Options:" << endl;
  cerr << "-dofs     : Write the transformations from each input surface to the procrustes mean." << endl;
  cerr << "-v        : verbose." << endl;
  cerr << "-m [file] : write mean shape to [file]" << endl;
  exit(0);
}

int main(int argc, char **argv)
{
  // This code does a Procustes alignment of a series of shapes using
  // scaling, rotations and translations output is all the aligned shapes

  int i, j, k, ok, writeDofs, writeMean, verbose;

  writeDofs = False;
  writeMean = False;
  verbose = False;

  irtkMatrix *matrices;
  irtkMatrix dummyMatrix(4, 4);
  irtkAffineTransformation *dof = new irtkAffineTransformation;
  double coord[3];
  char buf[300];
  char dofSuffix[] = "dof\0";

  if (argc < 3){
    usage();
  }
  
  // Number of datasets.
  int iNoOfDatasets = atoi(argv[1]);
  argv++;
  argc--;

  // Template for connectivity. Used when writing results.
  templateFileName = argv[1];

  int iNoOfLandmarks_1, iNoOfLandmarks_2;
  iNoOfLandmarks_1 = iNoOfLandmarks_2 = 0;

  irtkMatrix M_1;

  // Read in the names of the input files.
  inputFilenames = new char*[iNoOfDatasets];
  for (i = 0; i < iNoOfDatasets; i++){
    inputFilenames[i] = argv[1];
    argv++;
    argc--;
  }

  // Read in the names of the output files.
  outputFilenames = new char*[iNoOfDatasets];
  for (i = 0; i < iNoOfDatasets; i++){
    outputFilenames[i] = argv[1];
    argv++;
    argc--;
  }

  // Make space for transformation matrices.
  matrices = new irtkMatrix[iNoOfDatasets];
  for (i = 0; i < iNoOfDatasets; i++){
    matrices[i].Initialize(4, 4);
  }

  // Read optional arguments.
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofs") == 0)){
      argc--;
      argv++;
      writeDofs = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-v") == 0)){
      argc--;
      argv++;
      verbose = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-m") == 0)){
      argc--;
      argv++;
      writeMean = True;
      meanFilename = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-av") == 0)){
      argc--;
      argv++;
      alignVersion = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }

    if (ok == False){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  } 

 // Read all landmarks_1 in matrix M_1
  for (i = 0; i < iNoOfDatasets; i++){

    data1    = read(inputFilenames[i]);

    points1  = data1->GetPoints();
    vectors1 = data1->GetPointData()->GetVectors();

    if (i == 0){
      iNoOfLandmarks_1 = points1->GetNumberOfPoints();
      M_1.Initialize(3 * iNoOfLandmarks_1, iNoOfDatasets);
      cerr << "There are " << iNoOfDatasets
           << " datasets with "
           << points1->GetNumberOfPoints()
           << " landmarks." << endl;
    }

    if (verbose){
      cerr << "Including landmarks in " << inputFilenames[i] << endl;
    }

    for (j = 0; j < iNoOfLandmarks_1; j++){
      double p[3];
      
	points1->GetPoint(j, p);
	M_1(3 * j, i)   = p[0];
	M_1(3 * j + 1, i) = p[1];
	M_1(3 * j + 2, i) = p[2];
    }
  }

  // matrix of all aligned shapes
  irtkMatrix aligned_shapes(3 * iNoOfLandmarks_1, iNoOfDatasets); 
  
  // align each shape to first data set and calculate mean of aligned shapes
  irtkMatrix first_shape;

  first_shape.Initialize(3 * iNoOfLandmarks_1, 1);

  for (j = 0; j < 3 * iNoOfLandmarks_1; j++)
    first_shape(j, 0) = M_1(j, 0);
  
  // current mean pre alignment
  irtkMatrix current_mean;
  // initialize as first shape
  current_mean = first_shape;
  
  // shape we are currently aligning to aligned mean
  irtkMatrix current_shapetoalign;
  current_shapetoalign.Initialize(3 * iNoOfLandmarks_1, 1);
  
  // current shape after aligning to aligned mean
  irtkMatrix current_alignedshape;
  for (i = 0; i < iNoOfDatasets; i++)
  {
    // align ith shape
    for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
      current_shapetoalign(j, 0) = M_1(j, i);
    }

    //xxx
    current_alignedshape = align(first_shape, current_shapetoalign, iNoOfLandmarks_1, dummyMatrix);

    if (verbose){
      cerr << "aligned shape " << i << " to first shape" << endl;
    }
   
    // add to current mean
    for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
      current_mean(j, 0) += current_alignedshape(j, 0);
    }
  }

  // now have mean of aligned shapes
  current_mean *= (1.0 / iNoOfDatasets);
  
  float error2 = 0.0;

  for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
    error2 += (first_shape(j, 0) - current_mean(j, 0)) * (first_shape(j, 0) - current_mean(j, 0));
  }
  cerr << "error = " << error2/(3 * iNoOfLandmarks_1) << endl;; 
  
  // align mean to first data set

  //stores aligned mean that we use going into loop
  irtkMatrix previous_aligned_mean;

  //xxx
  previous_aligned_mean = align(first_shape, current_mean, iNoOfLandmarks_1, dummyMatrix);
  cerr << "aligned mean to first shape" << endl;
  
  irtkMatrix current_aligned_mean; //stores aligned mean calculated at end of loop
  int noits = 0;
  float error;
  do{
    // reinitialize current mean
    for (j = 0;j<3 * iNoOfLandmarks_1;j++){
      current_mean(j, 0) = 0.0;
    }
  
    // realign every shape with previous aligned mean 
    for (i = 0; i < iNoOfDatasets; i++){
    // align ith shape
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
        current_shapetoalign(j, 0) = M_1(j, i);
      }

      //xxx
      current_alignedshape = align(previous_aligned_mean, current_shapetoalign, iNoOfLandmarks_1, matrices[i]);

      for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
        aligned_shapes(j, i) = current_alignedshape(j, 0);
      }
      // add to current mean
      for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
        current_mean(j, 0)+=current_alignedshape(j, 0);
      }
    }
  
  
    current_mean *= (1.0 / iNoOfDatasets);

    //xxx
    current_aligned_mean = align(first_shape, current_mean, iNoOfLandmarks_1, dummyMatrix);


    for (i = 0; i < iNoOfDatasets; i++){
      matrices[i] = dummyMatrix * matrices[i];
    }

    error = 0;
    for (j = 0; j < 3 * iNoOfLandmarks_1; j++){
      error+= (current_aligned_mean(j, 0) - previous_aligned_mean(j, 0)) * (current_aligned_mean(j, 0) - previous_aligned_mean(j, 0));
    }

    cerr << "error = " << error << endl;
    noits++;

    previous_aligned_mean = current_aligned_mean;

    cerr << "no. of iterations = " << noits << endl;

  } while (error>0.0000001);

  
  // Write aligned shapes as vtk output
  // Read surface used as template.
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(templateFileName);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  for (j = 0; j < iNoOfDatasets; j++){

    for (i = 0; i < surface->GetNumberOfPoints(); i++){

      coord[0] = aligned_shapes(3 * i    , j);
      coord[1] = aligned_shapes(3 * i + 1, j);
      coord[2] = aligned_shapes(3 * i + 2, j);

      surface->GetPoints()->SetPoint(i, coord);
    }

    surface->Modified();
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(surface);
    if (verbose){
      cout << "Writing output " << j + 1 << " to " << outputFilenames[j] << endl;
    }

    writer->SetFileName(outputFilenames[j]);
    writer->Write();

    if (writeDofs == True){
      strcpy(buf, outputFilenames[j]);
      k = 0;
      while (buf[k] != '\0'){
        k++;
      }
      strcpy(&buf[k - 3], dofSuffix);
      if (verbose){
        cout << "Dof file is : " << buf << endl;
      }
      dof->PutMatrix(matrices[j]);
      dof->irtkTransformation::Write(buf);
    }

  }

  if (writeMean == True && meanFilename != NULL){

    for (i = 0; i < surface->GetNumberOfPoints(); i++){

      coord[0] = current_aligned_mean(3 * i    , 0);
      coord[1] = current_aligned_mean(3 * i + 1, 0);
      coord[2] = current_aligned_mean(3 * i + 2, 0);

      surface->GetPoints()->SetPoint(i, coord);
    }

    cout << "Writing mean to " << meanFilename << endl;

    surface->Modified();
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(surface);
    writer->SetFileName(meanFilename);
    writer->Write();
  }

  exit(0);

  ///////////////////// TESTING

  double ssq, meanPt[3], sd, maxSD;


  for (i = 0; i < surface->GetNumberOfPoints(); i++){

    coord[0] = coord[1] = coord[2] = 0;

    for (j = 0; j < iNoOfDatasets; j++){
      coord[0] += (M_1(3 * i    , j) / iNoOfDatasets);
      coord[1] += (M_1(3 * i + 1, j) / iNoOfDatasets);
      coord[2] += (M_1(3 * i + 2, j) / iNoOfDatasets);
    }

    surface->GetPoints()->SetPoint(i, coord);

  }

  surface->Modified();
  surface->Update();

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars->SetNumberOfComponents(1);

  maxSD = 0;


  for (i = 0; i < surface->GetNumberOfPoints(); i++){
    surface->GetPoint(i, meanPt);
    ssq = 0;

    for (j = 0; j < iNoOfDatasets; j++){
      ssq += (M_1(3 * i    , j) - meanPt[0]) * (M_1(3 * i    , j) - meanPt[0]);
      ssq += (M_1(3 * i + 1, j) - meanPt[1]) * (M_1(3 * i + 1, j) - meanPt[1]);
      ssq += (M_1(3 * i + 2, j) - meanPt[2]) * (M_1(3 * i + 2, j) - meanPt[2]);
    }

    sd = sqrt(ssq / iNoOfDatasets);
    scalars->InsertTuple1(i, sd);
    if (sd > maxSD){
      maxSD = sd;
    }

  }

  surface->GetPointData()->SetScalars(scalars);

  surface->Modified();
  surface->Update();

  cout << "Writing second mean with scalars to blah.vtk" << endl;
  cout << "Maximum SD = " << maxSD << endl;

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName("blah.vtk");
  writer->Write();


}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library"
       << endl;
}
#endif

