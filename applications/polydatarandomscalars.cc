#if (defined HAS_VTK)

#include <irtkImage.h>

//#include <nr.h>
//#include <time.h>
//#include <unistd.h>
#include <gsl/gsl_rng.h>

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
  cerr << " " << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  long i;

  int noOfPoints;
  bool ok;
  time_t seconds;
//  long ran2Seed;
//  long ran2initialSeed;

  if (argc < 3){
    usage();
  }

  // Prepare for random stuff.
//  pid_t procId = getpid();
//
//  seconds = time(NULL);
//  ran2Seed = seconds * procId;
//  if (ran2Seed < 1)
//  	ran2Seed = -1 * ran2Seed;
//
//  ran2initialSeed = -1 * ran2Seed;
//  cout << "Initial seed : " << ran2initialSeed << endl;
//
//  (void) ran2(&ran2initialSeed);

  gsl_rng * ranGen;
  const gsl_rng_type * ranGenType;
  gsl_rng_env_setup();
  ranGenType = gsl_rng_default;
  ranGen = gsl_rng_alloc (ranGenType);

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
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Input      : " << input_name << endl;
  cerr << "Output     : " << output_name << endl;

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData* input = reader->GetOutput();
  input->Update();

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
  writer->SetInput(input);
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
