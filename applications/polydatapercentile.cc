#include <irtkImage.h>
#include <nr.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#define MAXVALS 100

char *input_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << " polydatapercentile  [input] [N] [prctile_1] ... [prctile_N] <-options>" << endl;
  cerr << " Find percentiles of the scalars for a surface." << endl;
  cerr << " N : number of percentiles required (max " << MAXVALS << ")." << endl;
  cerr << " [prctile_1] ... [prctile_N] : percentiles required." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -q         : Quiet reporting; space separated percentiles only." << endl;
  cerr << " -name Name : Use scalars with given name." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  int i, numberOfPercentiles;
  int percentile[MAXVALS];
  int quiet = false;

  bool ok;
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
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
       argc--;
       argv++;
       scalar_name = argv[1];
       argc--;
       argv++;
       ok = true;
     }

    if (ok == false){
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

  vtkFloatArray *scalars;
  int ind;

  if (scalar_name == NULL){
    scalars = (vtkFloatArray*) input->GetPointData()->GetScalars();
  } else {
    scalars = (vtkFloatArray*) input->GetPointData()->GetArray(scalar_name, ind);

    if (ind == -1 || scalars == NULL){
      cerr << "Scalars unavailable with name " << scalar_name << endl;
      exit(0);
    }

  }

  // NR subroutines use 1-indexing.
  float *data = new float[1 + count];

  for (i = 0; i < count; ++i){
    data[1 + i] = scalars->GetTuple1(i);
  }

  sort(count, data);

  if (quiet == true){
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
