#ifdef HAS_VTK


#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>


char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: vtkscale [vtk] [scale factor] [vtkout]" << endl;
  exit(0);
}

int main(int argc, char **argv)
{

  int nx, ny, nz, n, i;
  double p[3], v[3], scale;

  if (argc < 4)
    usage();

  // Parse file names
  input_name  = argv[1];
  argc--;
  argv++;
  scale       = atof(argv[1]);
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  vtkStructuredGridReader *reader = vtkStructuredGridReader::New();
  reader->SetFileName(input_name);
  reader->Update();
  grid = reader->GetOutput();

  nx = grid->GetDimensions()[0];
  ny = grid->GetDimensions()[1];
  nz = grid->GetDimensions()[2];
  n = nx * ny * nz;

  // Get pointer to points and vectors
  vtkDataArray* vectors    = grid->GetPointData()->GetVectors();
  vtkDataArray* vectorsOut = vtkFloatArray::New();
  vectorsOut->SetNumberOfComponents(3);

  vtkPoints *points = vtkPoints::New();

  for (i = 0; i < n; ++i){
    p[0] = grid->GetPoint(i)[0];
    p[1] = grid->GetPoint(i)[1];
    p[2] = grid->GetPoint(i)[2];

    points->InsertNextPoint(p);

    v[0] = scale * vectors->GetTuple(i)[0];
    v[1] = scale * vectors->GetTuple(i)[1];
    v[2] = scale * vectors->GetTuple(i)[2];
    vectorsOut->InsertNextTuple(v);
  }

  // Allocate objects for vtkStructuredGrid format
  vtkStructuredGrid *gridOut = vtkStructuredGrid::New();
  vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();

//   cout << "Number of points  : " << points->GetNumberOfPoints(); 
//   cout << endl;
//   cout << "(nx, ny, nz) = " << "(" << nx << ", " << ny << ", " << nz << ")" << endl;
//   cout << "Number of vectors:  " << vectorsOut->GetNumberOfTuples() << endl;
 
  // Set structured grid
  gridOut->SetDimensions(nx, ny, nz);
  gridOut->SetPoints(points);
  gridOut->GetPointData()->SetVectors(vectorsOut);

//   gridOut->GetPoint(0, &p[0]);
//   cout << "gridOut->GetPoint(0, p[0], p[1], p[2]) " << p[0] << " " << p[1] << " " << p[2] << endl;

  // Write structured grid
  writer->SetInputData(gridOut);
  writer->SetFileName(output_name);
  writer->SetVectorsName("vectors");
  writer->Update();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
