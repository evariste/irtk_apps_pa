#include <irtkImage.h>
#include <nr.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#define MAXVALS 100

char *input_name = NULL;

void usage()
{
  cerr << " polydatapercentile  [input] [N] [prctile_1] ... [prctile_N] <-options>" << endl;
  cerr << " Find percentiles of the scalars for a surface." << endl;
  cerr << " N : number of percentiles required (max " << MAXVALS << ")." << endl;
  cerr << " [prctile_1] ... [prctile_N] : percentiles required." << endl;
  cerr << " Options:" << endl;
  cerr << " -q   : Quiet reporting; space separated percentiles only." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  int i, numberOfPercentiles;
  int percentile[MAXVALS];
  int quiet = False;

  int ok;
  int count, index;

  input_name  = argv[1];
  argc--;
  argv++;
  numberOfPercentiles = atoi(argv[1]);
  argc--;
  argv++;

  if (numberOfPercentiles > MAXVALS){
    cerr << "Cannot have more than " << MAXVALS << " percentiles!" << endl << endl;
    usage();
    exit(1);
  }

  for (i = 0; i < numberOfPercentiles; i++){
    percentile[i] = atoi(argv[1]);
    argc--;
    argv++;
  }

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = True;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  // Read surface
  vtkPolyData *input = vtkPolyData::New();
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();
  input = reader->GetOutput();

  // Find number of values.
  count = input->GetNumberOfPoints();

  if (count != input->GetPointData()->GetScalars()->GetNumberOfTuples()){
    cerr << "Size of scalars doesn't match number of points." << endl;
    exit(1);
  }

  // NR subroutines use 1-indexing.
  float *data = new float[1 + count];

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars = (vtkFloatArray*) input->GetPointData()->GetScalars();

  for (i = 0; i < count; ++i){
    data[1 + i] = scalars->GetTuple1(i);
  }

  sort(count, data);

  if (quiet == True){
    for (i = 0; i < numberOfPercentiles; i++){
      index = 1 + (int) round( (double) percentile[i] * (count - 1) / 100.0);
      cout << data[index];
      if (i < numberOfPercentiles-1)
        cout << " ";
    }
  } else {
    for (i = 0; i < numberOfPercentiles; i++){
      index = 1 + (int) round( (double) percentile[i] * (count - 1) / 100.0);
      cout << "  Percentile " << percentile[i] << " = " << data[index] << endl;
    }
  }

}
