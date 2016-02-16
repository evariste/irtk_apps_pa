#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;
char *template_name = NULL;

void usage()
{
  cerr << "Usage: matrix2polydata [irtkMatrix]  [vtkFileOut] <-template file>" << endl;
  cerr << "The points in [irtkMatrix] are written to vtk file [vtkFileOut]." << endl;
  cerr << "No connectivity is set unless a template vtk file is given from" << endl;
  cerr << "which the connectivity is copied." << endl;
  cerr << "The matrix data is expected to be in a single column containing" << endl;
  cerr << "3 * numberOfPoints entries (x_1 y_1 z_1 ... x_n y_n z_n)^T" << endl;
  cerr << "WARNING: If a template file is used, the normals may need recalculating. " << endl;
  cerr << "Clearly a template vtk file must have the correct number of points." << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, numberOfPoints;
  bool ok;
  int rows, cols;
  double p[3];
  irtkMatrix matrix;

  if (argc < 3){
    usage();
  }

  // Parse source and target point lists
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-template") == 0)){
      argc--;
      argv++;
      template_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  matrix.Read(input_name);

  rows = matrix.Rows();
  cols = matrix.Cols();

  if (cols != 1){
    cerr << "Matrix must have a single column." << endl;
    usage();
  }

  if (rows % 3 != 0){
    cerr << "Matrix must have a multiple of 3 rows." << endl;
    usage();
  }

  numberOfPoints = rows / 3;

  vtkPolyData *polys = vtkPolyData::New();
  vtkPoints   *points  = vtkPoints::New();

  if (template_name != NULL){

    // Read the input points.
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(template_name);
    reader->Update();
    polys = reader->GetOutput();

    if (polys->GetNumberOfPoints() != numberOfPoints) {
      cerr << "Template file must have same number of points as input matrix." << endl;
      usage();
    }

  }

  for (i = 0; i < numberOfPoints; i++){
    p[0] = matrix(3 * i    , 0);
    p[1] = matrix(3 * i + 1, 0);
    p[2] = matrix(3 * i + 2, 0);
    points->InsertPoint(i, p);
  }

  polys->SetPoints(points);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(polys);
  writer->SetFileName(output_name);
  writer->Write();

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
