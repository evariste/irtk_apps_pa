#if (defined HAS_VTK)

// Needs to have at least VTK 5.2

#include <irtkImage.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>

#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkDecimatePro.h>

char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: polydataeuclidean [input] [output] <options>\n" << endl;
  cerr << "Expected value of Euclidean distance between points on surface. Randomly sampled a certain number of times. (see polydatageodesic). " << endl;
  cerr << "Options:" << endl;
  cerr << "-reps [N]  default=100." << endl;
  exit(1);
}

double mean(vtkFloatArray *scalars){
  int n, i;
  double sum;

  n = scalars->GetNumberOfTuples();
  sum = 0.0;

  for (i = 0; i < n; ++i){
    sum += scalars->GetComponent(i,0);
  }

  return sum / ((double) n);
}

int main(int argc, char *argv[])
{
  bool ok;
  int i;
  int numberOfPoints;
  int id1, id2;

    int reps = 100;
  double p1[3], p2[3];

  if (argc < 2){
    usage();
  }

  // Parse source and target point lists
  in_name  = argv[1];
  argc--;
  argv++;


  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-XX") == 0)){
      argc--;
      argv++;
      // Do stuff and possible increment argv etc
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-reps") == 0)){
      argc--;
      argv++;
      reps = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
      argc--;
      argv++;
      out_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Prepare for random stuff.
  gsl_rng * ranGen;
  const gsl_rng_type * ranGenType;
  gsl_rng_env_setup();
  ranGenType = gsl_rng_default;

  ranGen = gsl_rng_alloc (gsl_rng_mt19937);

  timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long init = tv.tv_usec;
  gsl_rng_set(ranGen, init);


  // Read the input surface.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);
  reader->Update();


  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInputData(reader->GetOutput());
  // Setting tolerance to 0 means that vtkMergePoints is used instead.
  cleaner->SetTolerance(0.0);
  cleaner->PointMergingOn();
  cleaner->Update();

  vtkTriangleFilter *triFilter2 = vtkTriangleFilter::New();
  triFilter2->SetInputData(cleaner->GetOutput());
  triFilter2->Update();

  // Select a pair of vertices at random.
  numberOfPoints = triFilter2->GetOutput()->GetNumberOfPoints();
  cout << "Points in triangulated, cleaned input mesh: " << numberOfPoints << endl;

  vtkPolyData *temp = vtkPolyData::New();
  temp = triFilter2->GetOutput();

  double runningTotal = 0.0;

  for (i = 0; i < reps; ++i){
    id1 = (int) floor(numberOfPoints * gsl_rng_uniform(ranGen));

	  temp->GetPoint(id1, p1);

	  for (id2 = 0; id2 < numberOfPoints; ++id2){
	    temp->GetPoint(id2, p2);
	    runningTotal += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
	  }
  }

  cout << runningTotal / ((double) reps * numberOfPoints) << endl;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif



