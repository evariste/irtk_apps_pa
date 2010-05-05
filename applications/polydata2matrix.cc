#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: polydata2matrix [vtkFile] [irtkMatrixOut]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, numberOfPoints;
  double p[3];

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

  // Read the input points.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();

  numberOfPoints = input->GetNumberOfPoints();

  irtkMatrix m(3 * numberOfPoints, 1);

  for (i = 0; i < numberOfPoints; i++){
    input->GetPoint(i, p);
    m(3 * i    , 0) = p[0];
    m(3 * i + 1, 0) = p[1];
    m(3 * i + 2, 0) = p[2];
  }

  m.Write(output_name);

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
