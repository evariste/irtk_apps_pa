#if (defined HAS_VTK)

#include <irtkImage.h>

//#include <nr.h>
//#include <time.h>
//#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>

#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

char *input_name  = NULL;
char *output_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << " " << endl;
  cerr << " Usage: polydatarandomscalars [input] [output]  <options>" << endl;
  cerr << " " << endl;
  cerr << " Random scalars given to input surface.  Uniform distribution in [0,1]." << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -name name   : Give the output scalars the name 'name'" << endl;
  cerr << " -seed M      : Seed for random number generator." << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  long i;

  int noOfPoints;
  bool ok;

  if (argc < 3){
    usage();
  }

  // Prepare for random stuff.

  gsl_rng * ranGen;
//  const gsl_rng_type * ranGenType;
  gsl_rng_env_setup();
//  ranGenType = gsl_rng_default;
  ranGen = gsl_rng_alloc (gsl_rng_mt19937);

  timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long init = tv.tv_usec;

  // Parse arguments.
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
       argc--;
       argv++;
       scalar_name = argv[1];
       argc--;
       argv++;
       ok = true;
     }
    if ((ok == false) && (strcmp(argv[1], "-seed") == 0)) {
      argc--;
      argv++;
      init = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  gsl_rng_set(ranGen, init);

  cerr << "Input      : " << input_name << endl;
  cerr << "Output     : " << output_name << endl;

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData* input = reader->GetOutput();

  noOfPoints = input->GetNumberOfPoints();

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(noOfPoints);

  for (i = 0; i < noOfPoints; ++i){
//    scalars->SetTuple1(i, ran2(&ran2Seed));
    scalars->SetTuple1(i, gsl_rng_uniform(ranGen));
  }

  if (scalar_name == NULL){
 	  scalars->SetName("Random");
  } else {
  	scalars->SetName(scalar_name);
  }

  input->GetPointData()->AddArray(scalars);

  // Write the result.
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(input);
  writer->SetFileName(output_name);
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
