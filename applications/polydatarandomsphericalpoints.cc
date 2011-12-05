#ifdef HAS_VTK


#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include <nr.h>
#include <time.h>
#include <unistd.h>



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

  // Prepare for random stuff.
  pid_t procId = getpid();

  seconds = time(NULL);
  ran2Seed = seconds * procId;
  if (ran2Seed < 1)
  	ran2Seed = -1 * ran2Seed;

  ran2initialSeed = -1 * ran2Seed;
  (void) ran2(&ran2initialSeed);




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
  	theta = 2 * M_PI * ran2(&ran2Seed);
  	phi   = acos((2 * ran2(&ran2Seed)) -  1.0);
  	pt[0] = cos(theta) * sin(phi);
  	pt[1] = sin(theta) * sin(phi);
  	pt[2] = cos(phi);

  	randPoints->SetPoint(i, pt);
  }

  //u = rand(1,n);
  //v = rand(1,n);
  //
  //theta = 2 * pi * u;
  //phi   = acos((2 * v) - 1);
  //
  //pts = [ (cos(theta) .* sin(phi)) ; (sin(theta) .* sin(phi)) ; cos(phi) ];


  output->SetPoints(randPoints);
  output->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(output);
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


