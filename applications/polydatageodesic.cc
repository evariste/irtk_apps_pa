#if (defined HAS_VTK)

// Needs to have at least VTK 5.2

#include <irtkImage.h>

#include <gsl/gsl_rng.h>
#include <sys/time.h>


#include <vtkIdList.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkDecimatePro.h>



char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: polydatageodesic [input] <options>\n" << endl;
  cerr << "Expected geodesic distance over surface. Estimated over a number of random points sampled on surface (see polydataeuclidean)." << endl;
  cerr << "Options:" << endl;
  cerr << "   -output [filename]" << endl;
  cerr << "   -reps [N] (default 100)" << endl;
  cerr << "   -seed M   Seed for random number generator." << endl;
  exit(1);
}

double mean(vtkDoubleArray *scalars){

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
  int id1;

  int reps = 100;
  double *vals;

  if (argc < 2){
    usage();
  }

  timeval tv;
  gettimeofday(&tv, NULL);
  unsigned long init = tv.tv_usec;

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
    if ((ok == false) && (strcmp(argv[1], "-seed") == 0)) {
      argc--;
      argv++;
      init = atoi(argv[1]);
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
  gsl_rng_env_setup();

  ranGen = gsl_rng_alloc (gsl_rng_mt19937);
  gsl_rng_set(ranGen, init);



  // Read the input surface.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);
  reader->Update();

  vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
  triFilter->SetInputData(reader->GetOutput());

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInputData(reader->GetOutput());
  // Setting tolerance to 0 means that vtkMergePoints is used instead.
  cleaner->SetTolerance(0.0);//(0.005);
  cleaner->PointMergingOn();
  cleaner->Update();


  vtkTriangleFilter *triFilter2 = vtkTriangleFilter::New();
  triFilter2->SetInputData(cleaner->GetOutput());
  triFilter2->Update();

  vtkDijkstraGraphGeodesicPath *path = vtkDijkstraGraphGeodesicPath::New();
  path->SetInputData(triFilter2->GetOutput());

  numberOfPoints = triFilter2->GetOutput()->GetNumberOfPoints();

  cout << "Points in triangulated, cleaned input mesh: " << numberOfPoints << endl;

  vtkPolyData *temp = vtkPolyData::New();
  vtkDoubleArray *dists = vtkDoubleArray::New();
  dists->SetNumberOfComponents(1);
  dists->SetNumberOfTuples(numberOfPoints);

  vals = new double[reps];

  double runningTotal = 0.0;

  for (i = 0; i < reps; ++i){

    id1 = (int) floor(numberOfPoints * gsl_rng_uniform(ranGen));

    // Give the chosen vertex to the shortest path filter.
    path->SetStartVertex(id1);
    path->SetEndVertex(0);
    // Calculate geodesics for whole surface.
    path->StopWhenEndReachedOff();
    path->Update();

    path->GetCumulativeWeights(dists);

    runningTotal += mean(dists);
  }

  cout << runningTotal / ((double) reps) << endl;

  if (out_name != NULL){
    temp = cleaner->GetOutput();
    temp->GetPointData()->AddArray(dists);

    // Write the output with the scalar.
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(out_name);
    writer->SetInputData(temp);
    writer->Update();
    writer->Write();
  }

  delete [] vals;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif


