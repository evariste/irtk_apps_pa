#ifdef HAS_VTK

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include <gsl/gsl_rng.h>
#include <sys/time.h>




// Default filenames
char *output_name = NULL;

void usage()
{
  cerr << "Usage: " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i;
  bool ok;
  // default
  int nPts = 1000;



  time_t seconds;
  long ran2Seed;
  long ran2initialSeed;


  gsl_rng * r; 
  gsl_rng_env_setup();

  r = gsl_rng_alloc (gsl_rng_mt19937);

  timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long init = tv.tv_usec;
  gsl_rng_set(r, init);

  // Check command line
  if (argc < 2) {
    usage();
  }

  output_name = argv[1];
  argc--;
  argv++;

  vtkPolyData *output = vtkPolyData::New();


  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-n") == 0)) {
      argc--;
      argv++;
      nPts = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPoints *randPoints = vtkPoints::New();

  randPoints->SetNumberOfPoints(nPts);

  double theta;
  double phi;
  double pt[3];


  for (i = 0; i < nPts; i++) {


	theta = 2 * M_PI * gsl_rng_uniform(r);
  	phi   = acos((2 * gsl_rng_uniform(r)) -  1.0);

	pt[0] = cos(theta) * sin(phi);
  	pt[1] = sin(theta) * sin(phi);
  	pt[2] = cos(phi);

  	randPoints->SetPoint(i, pt);
  }

  output->SetPoints(randPoints);
  output->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(output);
  writer->SetFileName(output_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif


