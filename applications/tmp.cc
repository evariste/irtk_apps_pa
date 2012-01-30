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
#include <vtkDijkstraGraphGeodesicPath2.h>
#include <vtkDecimatePro.h>



char *in_name  = NULL, *out_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << "Usage: polydatageodesic [input] <options>\n" << endl;
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "-output [filename]" << endl;
  cerr << "-reps [N]" << endl;
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
	////////////////////////////////////////////// Log def stuff
	//	int i;
	//
	//	if (argc < 3){
	//		usage();
	//	}
	//
	//	// dof file with represents a velocity field.
	//  dofin_name  = argv[1];
	//  argc--;
	//  argv++;
	//  dofout_name  = argv[1];
	//  argc--;
	//  argv++;
	//
	//  irtkTransformation *transform = irtkTransformation::New(dofin_name);
	//
	//  irtkMultiLevelFreeFormTransformation *mffd =
	//  		new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *) transform));
	//
	//  irtkFreeFormTransformation3D *ffd =
	//  		dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));
	//
	//  irtkTransformation *transform2 = irtkTransformation::New(dofin_name);
	//
	//  irtkMultiLevelFreeFormTransformation *mffd2 =
	//  		new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *) transform2));
	//
	//  irtkFreeFormTransformation3D *ffd2 =
	//  		dynamic_cast<irtkFreeFormTransformation3D *>(mffd2->GetLocalTransformation(0));
	//
	//  ffd->Print();
	//  cout << ffd << endl << endl;
	//
	////  ffd2->Print();
	////  cout << ffd2 << endl << endl;
	//
	//  ffd->ExponentiateEuler();
	//
	//  for (i = 0; i < ffd2->NumberOfDOFs(); ++i)
	//  	ffd2->Put(i, -1 * ffd2->Get(i));
	//
	//
	////  ffd->Log();
	//  mffd->irtkTransformation::Write(dofout_name);
	//
	//  ffd2->ExponentiateEuler();
	//  mffd2->irtkTransformation::Write("inv.dof");
	//
	//
	//
	////  int a[3] = {23, 3, -2};
	////  for (i = 0; i < 3; ++i){
	////  	f(a[i]);
	////  }
	////
	////  for (i = 0; i < 3; ++i){
	////  	cout << a[i] << endl;
	////
	////  }



  bool ok;
  int i;
  int numberOfPoints;
  int id1, id2;
  time_t seconds;
  long ran2Seed;
  long ran2initialSeed;
  double val;
  double range[2];
  int reps = 100;
  double *vals;
  double alpha;

  id1 = -1;
  id2 = -1;
  alpha = 0.0;

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
    if ((ok == false) && (strcmp(argv[1], "-ids") == 0)){
      argc--;
      argv++;
      id1 = atoi(argv[1]);
      argc--;
      argv++;
      id2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-alpha") == 0)){
      argc--;
      argv++;
      alpha = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-scalars") == 0)){
      argc--;
      argv++;
      scalar_name = argv[1];
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
  cleaner->SetTolerance(0.0);//(0.005);
  cleaner->PointMergingOn();
  cleaner->Update();


  vtkTriangleFilter *triFilter2 = vtkTriangleFilter::New();
  triFilter2->SetInput(cleaner->GetOutput());
  triFilter2->Update();

  vtkPolyData *cleanedPolys = vtkPolyData::New();
  cleanedPolys = triFilter2->GetOutput();
  cleanedPolys->Update();


  if (scalar_name != NULL){
    vtkFloatArray *scalars;// = vtkFloatArray::New();
    int ind;
    scalars = (vtkFloatArray*) cleanedPolys->GetPointData()->GetArray(scalar_name, ind);
    cleanedPolys->GetPointData()->SetActiveScalars(scalar_name);
    cleanedPolys->Update();

    if (ind == -1 || scalars == NULL){
      cerr << "Scalars unavailable with name " << scalar_name << endl;
      exit(1);
    }
  }

  vtkDijkstraGraphGeodesicPath2 *path = vtkDijkstraGraphGeodesicPath2::New();
  path->SetInput(cleanedPolys);

  if (scalar_name != NULL){
  	path->SetEdgeWeightsFromScalars(1);
  }

  path->SetAlpha(alpha);

  vtkIndent indent = vtkIndent(1);



  numberOfPoints = cleaner->GetOutput()->GetNumberOfPoints();

  vtkPolyData *temp = vtkPolyData::New();
  vtkDoubleArray *dists = vtkDoubleArray::New();
  dists->SetNumberOfComponents(1);
  dists->SetNumberOfTuples(numberOfPoints);

  vals = new double[reps];

  double runningTotal = 0.0;

  if (id1 > -1 && id2 > -1 and id1 < numberOfPoints and id2 < numberOfPoints){
    path->SetStartVertex(id1);
    path->SetEndVertex(id2);
    // Calculate geodesics for whole surface.
    path->StopWhenEndReachedOn();

    path->PrintSelf(cout, indent);

    path->Update();

  } else {

  	for (i = 0; i < reps; ++i){

  		id1 = (int) floor(numberOfPoints * ran2(&ran2Seed));

  		// Give the chosen vertex to the shortest path filter.
  		path->SetStartVertex(id1);
  		path->SetEndVertex(0);
  		// Calculate geodesics for whole surface.
  		path->StopWhenEndReachedOff();
  		path->Update();

  		// Assign the scalars (distances) obtained from the shortest path filter.
  		//    dists = (vtkFloatArray*)path->GetGeodesicLength();

  		//    path->GetAllLengths(dists);
  		path->GetCumulativeWeights(dists);

  		//    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
  		//    pd = path->GetOutput();
  		//    pd->Update();
  		//
  		//    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
  		//    wr->SetInput(pd);
  		//    wr->SetFileName("bla.vtk");
  		//    wr->Update();
  		//    wr->Write();
  		//    exit(0);
  		//


  		runningTotal += mean(dists);
  	}

  	cout << runningTotal / ((double) reps) << endl;

  }

  if (out_name != NULL){
//    temp = cleaner->GetOutput();
//    temp->GetPointData()->AddArray(dists);
//    temp->Update();

    // Write the output with the scalar.
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(out_name);
    // writer->SetInput(temp);
    writer->SetInput(path->GetOutput());
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



