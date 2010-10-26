#if (defined HAS_VTK)

// Needs to have at least VTK 5.2

#include <irtkImage.h>
#include <nr.h>
#include <time.h>

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
  cerr << "" << endl;
  cerr << "Options:" << endl;
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
  bool ok;
  int i;
  int numberOfPoints;
  int id1, id2;
  time_t seconds;
  long ran2Seed;

  long ran2initialSeed;
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
  cleaner->SetTolerance(0.0);
  cleaner->PointMergingOn();
  cleaner->Update();

  vtkTriangleFilter *triFilter2 = vtkTriangleFilter::New();
  triFilter2->SetInput(cleaner->GetOutput());
  triFilter2->Update();

  // Select a pair of vertices at random.
  numberOfPoints = cleaner->GetOutput()->GetNumberOfPoints();

  vtkPolyData *temp = vtkPolyData::New();
  temp = cleaner->GetOutput();

  double runningTotal = 0.0;

  for (i = 0; i < reps; ++i){
	  id1 = (int) floor(numberOfPoints * ran2(&ran2Seed));

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



