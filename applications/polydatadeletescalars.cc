#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydatadeletescalars [input] [output] [scalarName]" << endl;
  cerr << "" << endl;
  exit(1);
}

char *input_name = NULL;
char *output_name = NULL;
char *scalar_name = NULL;

int main(int argc, char **argv)
{
  if (argc < 4){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;
  scalar_name  = argv[1];
  argc--;
  argv++;

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();

  vtkPolyData *surface = surface_reader->GetOutput();

  int i, noOfArrays;
  int success = False;

  noOfArrays = surface->GetPointData()->GetNumberOfArrays();

  for (i = 0; i < noOfArrays; i++){
    vtkFloatArray *scalars;

    scalars = (vtkFloatArray*) surface->GetPointData()->GetArray(i);

    if ( scalars->GetNumberOfComponents() != 1 )
      continue;

    if (strcmp(scalars->GetName(), scalar_name) == 0){
      cout << "Deleting scalars : " << scalars->GetName() << endl;
      surface->GetPointData()->RemoveArray(scalar_name);
      surface->Update();
      success = True;
      break;
    }
  }

  if (success == True){
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(surface);
    writer->SetFileName(output_name);
    writer->Write();
  } else {
    cerr << "No such scalars : " << scalar_name << endl;
  }



}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
