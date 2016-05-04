#if (defined HAS_VTK)


#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkHull.h>

char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: polydatahull [input] [output] <options>\n" << endl;
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "-levels n  : Number of levels to subdivide planes in hull." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int levels = 1;

  if (argc < 3){
    usage();
  }

  // Parse source and target point lists
  in_name  = argv[1];
  argc--;
  argv++;
  out_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-levels") == 0)){
      argc--;
      argv++;
      levels = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // Read the input surface.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();

  vtkHull *hull = vtkHull::New();
  hull->SetInputData(input);
  hull->AddRecursiveSpherePlanes(levels);

  hull->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(out_name);
  writer->SetInputData(hull->GetOutput());

  writer->Update();
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
