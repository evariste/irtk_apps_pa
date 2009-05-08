#if (defined HAS_VTK)


#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: polydatarecalculatenormals [input] [output] <-splittingOn>\n" << endl;
  cerr << "Default splitting is off." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  int splitting = False;

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
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-splittingOn") == 0)){
      argc--;
      argv++;
      splitting = True;
      ok = True;
    }
    if (ok == False){
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

  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();

  normalsFilter->SetInput(input);

  if (splitting == True){
    cerr << "Setting splitting to on" << endl;
    normalsFilter->SplittingOn();
  } else {
    cerr << "Setting splitting to off" << endl;
    normalsFilter->SplittingOff();
  }
  normalsFilter->Modified();
  normalsFilter->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(out_name);
  writer->SetInput(normalsFilter->GetOutput());
  writer->Modified();
  writer->SetFileTypeToBinary();

  writer->Update();
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
