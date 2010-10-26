#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkPolyData.h>

#include <vtkPolyDataAlgorithm.h>

#include <vtkPolyDataNormals.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydatareflect [in] [out] [-x|-y|-z]" << endl;
  exit(1);
}

char *input_name = NULL, *output_name = NULL, *dof_name  = NULL;
int main(int argc, char **argv)
{

  if (argc < 3){
    usage();
  }

  int i;
  bool ok;

  int reflectX = false;
  int reflectY = false;
  int reflectZ = false;
  double coord[3];

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      reflectX = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      reflectY = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-z") == 0)){
      argc--;
      argv++;
      reflectZ = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  for (i = 0; i < surface->GetNumberOfPoints(); i++){

    (surface->GetPoints())->GetPoint(i, coord);

    if (reflectX == true)
      coord[0] = -1.0 * coord[0];

    if (reflectY == true)
      coord[1] = -1.0 * coord[1];

    if (reflectZ == true)
      coord[2] = -1.0 * coord[2];

    surface->GetPoints()->SetPoint(i, coord);
  }

  surface->Modified();

  //Update the normals to reflect the new points.
  vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
  normals->SplittingOff();
  normals->SetInput(surface);
  normals->Update();
  surface = normals->GetOutput();
  surface->Modified();


  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(output_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
