#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydataaddtoscalars [input] [value] [output]"<<endl;
//   cerr << "Options:" << endl;
  cerr << "Add [value] to scalars in [input] polydata." << endl;
  cerr << "Write to [output]" << endl;
  exit(1);
}

int main(int argc, char **argv ){

  if (argc < 4 ){
    usage();
  }

  char *input_name;
  char *output_name;
  int i;
  double value, val;
  int noOfPoints;

  input_name = argv[1];
  argv++;
  argc--;

  value = atof(argv[1]);
  argv++;
  argc--;

  output_name = argv[1];
  argv++;
  argc--;


  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  vtkPolyData *polys = vtkPolyData::New();

  reader->SetFileName(input_name);
  reader->Modified();
  reader->Update();

  polys= reader->GetOutput();

  noOfPoints = polys->GetNumberOfPoints();
  cout << "No of points " << noOfPoints << endl;

  vtkFloatArray *scalars = vtkFloatArray::New();

  scalars = (vtkFloatArray*) polys->GetPointData()->GetScalars();

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    val += value;
    scalars->SetTuple1(i, val);
  }

  scalars->Modified();

  polys->GetPointData()->SetScalars(scalars);
  polys->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(output_name);
  writer->SetInput(polys);
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
