#if (defined HAS_VTK)

// Needs to have at least VTK 5.2

#include <irtkImage.h>

#include <nr.h>
#include <time.h>


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
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "-output [filename]" << endl;
  cerr << "-reps [N]" << endl;
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
  int ok, i;
  int numberOfPoints;
  int id1, id2;
  time_t seconds;
  long ran2Seed;
  long ran2initialSeed;
  double val;
  double range[2];
  int reps = 100;
  double *vals;

  if (argc < 2){
    usage();
  }

  // Parse source and target point lists
  in_name  = argv[1];
  argc--;
  argv++;


  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-XX") == 0)){
      argc--;
      argv++;
      // Do stuff and possible increment argv etc
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-reps") == 0)){
      argc--;
      argv++;
      reps = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-output") == 0)){
      argc--;
      argv++;
      out_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }

    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Prepare for random stuff.
  seconds = time(NULL);
  ran2Seed = seconds;
  ran2initialSeed = -1 * ran2Seed;
  (void) ran2(&ran2initialSeed);

  // Read the input surface.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);

  vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
  triFilter->SetInput(reader->GetOutput());

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInput(reader->GetOutput());
  // Setting tolerance to 0 means that vtkMergePoints is used instead.
  cleaner->SetTolerance(0.0);//(0.005);
  cleaner->PointMergingOn();
  cleaner->Update();


  vtkTriangleFilter *triFilter2 = vtkTriangleFilter::New();
  triFilter2->SetInput(cleaner->GetOutput());
  triFilter2->Update();

  vtkDijkstraGraphGeodesicPath *path = vtkDijkstraGraphGeodesicPath::New();
  path->SetInput(triFilter2->GetOutput());

  numberOfPoints = cleaner->GetOutput()->GetNumberOfPoints();

  vtkPolyData *temp = vtkPolyData::New();
  vtkFloatArray *dists = vtkFloatArray::New();

  vals = new double[reps];

  double runningTotal = 0.0;

  for (i = 0; i < reps; ++i){

    id1 = (int) floor(numberOfPoints * ran2(&ran2Seed));

    // Give the chosen vertex to the shortest path filter.
    path->SetStartVertex(id1);
    path->SetEndVertex(0);
    // Calculate geodesics for whole surface.
    path->SetStopWhenEndReached(0);
    path->Update();

    // Assign the scalars (distances) obtained from the shortest path filter.
    dists = (vtkFloatArray*)path->Getd();

    runningTotal += mean(dists);
  }

  cout << runningTotal / ((double) reps) << endl;

  if (out_name != NULL){
    temp = cleaner->GetOutput();
    temp->GetPointData()->SetScalars(dists);
    temp->Update();

    // Write the output with the scalar.
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(out_name);
    writer->SetInput(temp);
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


