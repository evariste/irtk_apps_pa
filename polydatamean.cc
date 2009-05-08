#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>


char **input_names = NULL, *output_name = NULL;
vtkPolyData **inputs;
vtkPolyDataReader **readers;

void usage()
{
  cerr << "Usage: polydatamean [outFile] N [polys_1] ... [polys_N]" << endl;
  cerr << "Note that output mean normals (if present) may need recalculating." << endl;
  exit(1);
}

int main(int argc, char **argv)
{

  int number, i, j, pointCount;
  double p_in[3];
  double p_mean[3];

  if (argc < 4){
    usage();
  }

  output_name = argv[1];
  argv++;
  argc--;
  number = atoi(argv[1]);
  argv++;
  argc--;

  // Error checking.
  if (argc - 1 != number){
    cerr << "Number of data sets ("<< argc - 1 << ") is not equal to paramater given (" << number << ")" << endl;
    exit(1);
  }
  if (number < 1){
    cerr << "Please give at least one dataset." << endl;
    exit(1);
  }

  // Read data.
  cout << "Reading " << number << " datasets ..." << endl;
  input_names = new char* [number];
  for (i = 0; i < number; ++i){
    input_names[i] = argv[1];
    argc--;
    argv++;
  }

  inputs  = new vtkPolyData *[number];
  readers = new vtkPolyDataReader *[number];

  for (i = 0; i < number; i++){
    inputs[i] = vtkPolyData::New();
    readers[i] = vtkPolyDataReader::New();

    readers[i]->SetFileName(input_names[i]);
    readers[i]->Modified();
    readers[i]->Update();
    inputs[i] = readers[i]->GetOutput();
  }
  cout << " done." << endl;

  // Error checking.
  pointCount = inputs[0]->GetNumberOfPoints();
  cout << "Number of points : " << pointCount << endl;

  for (i = 0; i < number; i++){
    if (inputs[i]->GetNumberOfPoints() != pointCount){
      cerr << "All datasets must have the same number of points." << endl;
      exit(1);
    }
  }

  // Copy first dataset into mean.
  vtkPolyData *mean = vtkPolyData::New();
  mean = readers[0]->GetOutput();

  // Averaging.
  for (j = 0; j < pointCount; j++){
    p_mean[0] = p_mean[1] = p_mean[2] = 0;

    for (i = 0; i < number; i++){
      inputs[i]->GetPoint(j, p_in);
      p_mean[0] += p_in[0];
      p_mean[1] += p_in[1];
      p_mean[2] += p_in[2];
    }

    p_mean[0] /= number;
    p_mean[1] /= number;
    p_mean[2] /= number;

    mean->GetPoints()->SetPoint(j, p_mean);
  }

  mean->Modified();

  // Write output
  cout << "Writing output to " << output_name << endl;
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(mean);
  writer->SetFileName(output_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
