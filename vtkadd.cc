#ifdef HAS_VTK



#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>


char *input1_name = NULL, *input2_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: vtkadd [vtk1] [vtk2] [vtkout]" << endl;
  exit(0);
}

int main(int argc, char **argv)
{

  int nx, ny, nz, n, i;
  double x1, x2, y1, y2, z1, z2, p[3], v[3];
  double tol = 0.00001;

  if (argc < 4)
    usage();

  // Parse file names
  input1_name  = argv[1];
  argc--;
  argv++;
  input2_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  vtkStructuredGrid *grid1 = vtkStructuredGrid::New();
  vtkStructuredGridReader *reader1 = vtkStructuredGridReader::New();
  reader1->SetFileName(input1_name);
  reader1->Update();
  grid1 = reader1->GetOutput();
  vtkStructuredGrid *grid2 = vtkStructuredGrid::New();
  vtkStructuredGridReader *reader2 = vtkStructuredGridReader::New();
  reader2->SetFileName(input2_name);
  reader2->Update();
  grid2 = reader2->GetOutput();

  nx = grid1->GetDimensions()[0];
  ny = grid1->GetDimensions()[1];
  nz = grid1->GetDimensions()[2];
  n = nx * ny * nz;

  if (nx != grid2->GetDimensions()[0] || ny != grid2->GetDimensions()[1] || nz != grid2->GetDimensions()[2]){
    cerr << "Input grids must have the same dimensions." << endl;
    usage();
  }

  x1 = grid1->GetPoint(0)[0];
  y1 = grid1->GetPoint(0)[1];
  z1 = grid1->GetPoint(0)[2];
  x2 = grid1->GetPoint(n-1)[0];
  y2 = grid1->GetPoint(n-1)[1];
  z2 = grid1->GetPoint(n-1)[2];

  if (fabs(x1 - grid2->GetPoint(0)[0]) > tol   ||
      fabs(y1 - grid2->GetPoint(0)[1]) > tol   ||
      fabs(z1 - grid2->GetPoint(0)[2]) > tol   ||
      fabs(x2 - grid2->GetPoint(n-1)[0]) > tol ||
      fabs(y2 - grid2->GetPoint(n-1)[1]) > tol ||
      fabs(z2 - grid2->GetPoint(n-1)[2]) > tol){
    cerr << "Grid corners must be the same." << endl;
    usage();
  }

  // Get pointer to points and vectors
  vtkDataArray* vectors1   = grid1->GetPointData()->GetVectors();
  vtkDataArray* vectors2   = grid2->GetPointData()->GetVectors();
  vtkDataArray* vectorsOut = vtkFloatArray::New();
  vectorsOut->SetNumberOfComponents(3);

  vtkPoints *points = vtkPoints::New();

  for (i = 0; i < n; ++i){
    p[0] = grid1->GetPoint(i)[0];
    p[1] = grid1->GetPoint(i)[1];
    p[2] = grid1->GetPoint(i)[2];

    points->InsertNextPoint(p);

    v[0] = vectors1->GetTuple(i)[0] + vectors2->GetTuple(i)[0];
    v[1] = vectors1->GetTuple(i)[1] + vectors2->GetTuple(i)[1];
    v[2] = vectors1->GetTuple(i)[2] + vectors2->GetTuple(i)[2];
    vectorsOut->InsertNextTuple(v);
  }

  // Allocate objects for vtkStructuredGrid format
  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();

//   cout << "Number of points  : " << points->GetNumberOfPoints(); 
//   cout << endl;
//   cout << "(nx, ny, nz) = " << "(" << nx << ", " << ny << ", " << nz << ")" << endl;
//   cout << "Number of vectors:  " << vectorsOut->GetNumberOfTuples() << endl;
 
  // Set structured grid
  grid->SetDimensions(nx, ny, nz);
  grid->SetPoints(points);
  grid->GetPointData()->SetVectors(vectorsOut);

  grid->GetPoint(0, &p[0]);
//   cout << "grid->GetPoint(0, p[0], p[1], p[2]) " << p[0] << " " << p[1] << " " << p[2] << endl;

  // Write structured grid
  writer->SetInput(grid);
  writer->SetFileName(output_name);
  //  writer->SetFileTypeToBinary();
  writer->SetVectorsName("vectors");
  writer->Update();

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
