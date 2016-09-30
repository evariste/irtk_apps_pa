#ifdef HAS_VTK

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

#include <gsl/gsl_rng.h>
#include <sys/time.h>

// Default filenames
char *output_name = NULL;

void usage(std::string name)
{
  cerr << endl;
  cerr << "Usage: " << endl;
  cerr << "  " << name << "  [output] <options>" << endl;
  cerr << "   -n N    : Number of points (default = 1000)." << endl;
  cerr << "   -seed M : Seed for random number generator." << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i;
  bool ok;
  // default
  int nPts = 1000;

  std::string exeName;

  exeName = argv[0];

  gsl_rng * rng;
  gsl_rng_env_setup();

  rng = gsl_rng_alloc (gsl_rng_mt19937);

  timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long init = tv.tv_usec;

  // Check command line
  if (argc < 2) {
    usage(exeName);
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
    if ((ok == false) && (strcmp(argv[1], "-seed") == 0)) {
      argc--;
      argv++;
      init = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage(exeName);
    }
  }

  // Seed the random number generator.
  gsl_rng_set(rng, init);

  vtkPoints *randPoints = vtkPoints::New();

  randPoints->SetNumberOfPoints(nPts);

  double theta;
  double phi;
  double pt[3];

  for (i = 0; i < nPts; i++) {
    theta = 2 * M_PI * gsl_rng_uniform(rng);
    phi   = acos((2 * gsl_rng_uniform(rng)) -  1.0);

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


