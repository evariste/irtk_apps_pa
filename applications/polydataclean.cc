#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>

void usage(){
  cerr << "polydataclean [input] [output] <-options>"<<endl;
  cerr << "Options:" << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv ){

  if (argc <3 ){
    usage();
  }

  bool ok;
  char *input_name;
  char *output_name;

  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-option") == 0)){
      argc--;
      argv++;
//stuff
      ok  = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  vtkPolyData* inputPoly;
  vtkPolyDataReader* reader;
  vtkCleanPolyData* cleaner;

  inputPoly = vtkPolyData::New();
  reader    = vtkPolyDataReader::New();
  cleaner   = vtkCleanPolyData::New();

  reader->SetFileName(input_name);
  reader->Modified();
  reader->Update();

  inputPoly = reader->GetOutput();

  cleaner->SetInput(inputPoly);
  cleaner->Modified();
  cleaner->Update();

  // Write out.
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(output_name);
  writer->SetInput(cleaner->GetOutput());
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
