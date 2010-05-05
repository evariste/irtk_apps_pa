#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkPolyDataReader.h>

void usage(){
  cerr << "polydatalistscalars [polydata]" << endl;
  cerr << "" << endl;
  exit(1);
}

char *input_name = NULL;

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();

  vtkPolyData *surface = surface_reader->GetOutput();

  int i, noOfArrays;
  noOfArrays = surface->GetPointData()->GetNumberOfArrays();

  for (i = 0; i < noOfArrays; i++){
    vtkFloatArray *scalars;

    scalars = (vtkFloatArray*) surface->GetPointData()->GetArray(i);

    if ( scalars->GetNumberOfComponents() != 1 )
      continue;

    cout << scalars->GetName() << endl;
  }


}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
