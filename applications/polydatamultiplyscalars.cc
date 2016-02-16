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
#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>

char *scalar_name = NULL;
char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << " polydatamultiplyscalars [input] [factor] [output] <-options>"<<endl;
//   cerr << "Options:" << endl;
  cerr << " Multiply the scalars in the [input] polydata by [factor]." << endl;
  cerr << " Write to [output]" << endl;
  cerr << " " << endl;
  cerr << " Options: " << endl;
  cerr << " " << endl;
  cerr << " -name ScalarsName : Multiply the scalars with the given name." << endl;
  
  exit(1);
}

int main(int argc, char **argv ){

  if (argc < 4 ){
    usage();
  }

  int i;
  bool ok;
  double factor, val;
  int noOfPoints;

  input_name = argv[1];
  argv++;
  argc--;

  factor = atof(argv[1]);
  argv++;
  argc--;

  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)) {
      argc--;
      argv++;
      scalar_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  vtkPolyData *polys = vtkPolyData::New();

  reader->SetFileName(input_name);
  reader->Modified();
  reader->Update();

  polys= reader->GetOutput();

  noOfPoints = polys->GetNumberOfPoints();
  cout << "No of points " << noOfPoints << endl;

  vtkFloatArray *scalars;
  if (scalar_name != NULL){
    scalars = (vtkFloatArray*) polys->GetPointData()->GetArray(scalar_name);
    if (scalars == NULL){
      cerr << "Cannot retrieve scalars : " << scalar_name << endl;
      exit(1);
    }
  } else {
    scalars = (vtkFloatArray*) polys->GetPointData()->GetScalars();
    cerr << "Using scalars :  " << scalars->GetName() << endl;
    scalar_name = scalars->GetName();
  }


  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    val *= factor;
    scalars->SetTuple1(i, val);
  }

  scalars->Modified();

  polys->GetPointData()->AddArray(scalars);
  polys->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(output_name);
  writer->SetInputData(polys);
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
