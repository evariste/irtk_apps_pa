#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkFloatArray.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointLocator.h>

#include <vtkTriangleFilter.h>
#include <vtkMath.h>

char *input_nameA = NULL;
char *input_nameB = NULL;

void usage()
{
  cerr << " polydatasurfacedistance [inputA] [inputB]" << endl;
  cerr << " " << endl;
  cerr << " Estimate the 'currents' type distance between a pair" << endl;
  cerr << " of surfaces. See the Glaunes IPMI 2009 paper and " << endl;
  cerr << " Linh Ha, MICCAI, 2010." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
	int i, j, ok;
  double searchRadius;
  double a[3], b[3], c[3];
  double e1[3], e2[3], e3[3];
  double n1[3], n2[3];
  double volA, volB, rA, rB, sigmaKer, kernelFraction;
  int nearestPtCount;
  int noOfFacesA, noOfFacesB;
  double boundsA[6], boundsB[6];
  double constA, constB;

  kernelFraction = 0.05;


	if (argc < 2){
    usage();
  }

  input_nameA = argv[1];
  argc--;
  argv++;
  input_nameB = argv[1];
  argc--;
  argv++;


  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((!ok) && (strcmp(argv[1], "-bla") == 0)) {
    	//      argc--;
    	//      argv++;
    	//      stuff
    	// maybe
    	//      argc--;
    	//      argv++;
    	//      ok = True;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }


  // Read surface
  vtkPolyDataReader *surface_readerA = vtkPolyDataReader::New();
  surface_readerA->SetFileName(input_nameA);
  surface_readerA->Modified();
  surface_readerA->Update();

  vtkPolyDataReader *surface_readerB = vtkPolyDataReader::New();
  surface_readerB->SetFileName(input_nameB);
  surface_readerB->Modified();
  surface_readerB->Update();


  vtkPolyData *readerAout;
  readerAout = surface_readerA->GetOutput();
  readerAout->Update();
  readerAout->BuildCells();
  readerAout->BuildLinks();

  vtkPolyData *readerBout;
  readerBout = surface_readerB->GetOutput();
  readerBout->Update();
  readerBout->BuildCells();
  readerBout->BuildLinks();


  // Derive a new polydata set with points equal to the
  // centers of the faces of the input surfaces and vectors
  // associated with each point that are the normals of the faces.

  vtkTriangleFilter *triFilterA = vtkTriangleFilter::New();
  triFilterA->SetInput(readerAout);
  triFilterA->Update();
  vtkPolyData *triAout = vtkPolyData::New();
  triAout = triFilterA->GetOutput();
  triAout->BuildCells();
  triAout->BuildLinks();

  vtkTriangleFilter *triFilterB = vtkTriangleFilter::New();
  triFilterB->SetInput(readerBout);
  triFilterB->Update();
  vtkPolyData *triBout = vtkPolyData::New();
  triBout = triFilterB->GetOutput();
  triBout->BuildCells();
  triBout->BuildLinks();


  triAout->GetBounds(boundsA);
  triBout->GetBounds(boundsB);

  if (boundsA[0] > boundsB[1] ||
  		boundsB[0] > boundsA[1] ||
  		boundsA[2] > boundsB[3] ||
  		boundsB[2] > boundsA[3] ||
  		boundsA[4] > boundsB[5] ||
  		boundsB[4] > boundsA[5]){
  	cout << "Surfaces' bounding boxes do not intersect." << endl;
  	exit(1);
  }



  vtkCellArray *facesA = triAout->GetPolys();
  noOfFacesA = facesA->GetNumberOfCells();

  vtkCellArray *facesB = triBout->GetPolys();
  noOfFacesB = facesB->GetNumberOfCells();


  vtkPoints *centres = vtkPoints::New();
  centres->SetNumberOfPoints(noOfFacesA + noOfFacesB);

  vtkFloatArray *normals = vtkFloatArray::New();
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(noOfFacesA + noOfFacesB);

  vtkIdType nptsForFace = 0;
  vtkIdType *ptIdsForFace;


  facesA->InitTraversal();
  for (i = 0; i < noOfFacesA; ++i){

  	facesA->GetNextCell(nptsForFace, ptIdsForFace);

  	if (nptsForFace != 3){
  		cerr << "Error: faces must have three points each." << endl;
  		exit(0);
  	}

  	triAout->GetPoint(ptIdsForFace[0], a);
  	triAout->GetPoint(ptIdsForFace[1], b);
  	triAout->GetPoint(ptIdsForFace[2], c);

  	for (j = 0; j < 3; ++j){
  		e1[j] = b[j] - a[j];
  		e2[j] = c[j] - a[j];
  	}

  	vtkMath::Cross(e1, e2, e3);
  	normals->SetTuple3(i, e3[0], e3[1], e3[2]);

  	// Get centroid.
  	for (j = 0; j < 3; ++j){
  		a[j] += (b[j] + c[j]);
  		a[j] /= 3.0;
  	}

  	centres->SetPoint(i, a);
  }

  facesB->InitTraversal();
  for (i = noOfFacesA; i < noOfFacesA + noOfFacesB; ++i){

  	facesB->GetNextCell(nptsForFace, ptIdsForFace);

  	if (nptsForFace != 3){
  	  		cerr << "Error: faces must have three points each." << endl;
  	  		exit(0);
  	}

  	triBout->GetPoint(ptIdsForFace[0], a);
  	triBout->GetPoint(ptIdsForFace[1], b);
  	triBout->GetPoint(ptIdsForFace[2], c);

  	for(j = 0; j < 3; ++j){
  		e1[j] = b[j] - a[j];
  		e2[j] = c[j] - a[j];
  	}

  	// Cross product is skew symmetric for the second surface, we do e2 x e1 while we did
  	// e1 x e2 for the first surface - i.e.  the sense is negated.
  	vtkMath::Cross(e2, e1, e3);
  	normals->SetTuple3(i, e3[0], e3[1], e3[2]);

  	// Get centroid.
  	for (j = 0; j < 3; ++j){
  		a[j] += (b[j] + c[j]);
  		a[j] /= 3.0;
  	}

  	centres->SetPoint(i, a);
  }

  normals->SetName("faceNormals");

  vtkPolyData *combined = vtkPolyData::New();

  combined->SetPoints(centres);
  combined->GetPointData()->AddArray(normals);
  combined->Update();

//  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//  writer->SetInput(combined);
//  writer->SetFileTypeToBinary();
//  writer->SetFileName("temp.vtk");
//  writer->Write();
//  exit(0);

  vtkPointLocator *point_locator = vtkPointLocator::New();
  point_locator->SetNumberOfPointsPerBucket(5);
  point_locator->SetDataSet(combined);
  point_locator->BuildLocator();


  // Need to establish a suitable radius for the kernel on centre
  // to centre distances.

  volA = (boundsA[1] - boundsA[0]) *
				 (boundsA[3] - boundsA[2]) *
  		   (boundsA[5] - boundsA[4]);

  volB = (boundsB[1] - boundsB[0]) *
  		   (boundsB[3] - boundsB[2]) *
  		   (boundsB[5] - boundsB[4]);

  // Assume points distributed on spheres.
  rA = pow(volA, (1.0/3.0)) * 3.0 / 4.0 / M_PI;
  rB = pow(volB, (1.0/3.0)) * 3.0 / 4.0 / M_PI;

  sigmaKer = 0.5 * (rA + rB) * kernelFraction;


  // Constants for kernel
  constA = 1 / sigmaKer / sqrt(2 * M_PI);
  constB = -1.0 / sigmaKer / sigmaKer;

  searchRadius = 2.5 * sigmaKer;


  vtkIdList *nearestPtIDs = vtkIdList::New();

  double val, pt2ptDistSq, totalDist;

  // The distance we seek to measure and return.
  totalDist = 0.0;

  for (i = 0; i < noOfFacesA + noOfFacesB; ++i){
  	combined->GetPoint(i, a);
  	point_locator->FindPointsWithinRadius(searchRadius, a, nearestPtIDs);

  	nearestPtCount = nearestPtIDs->GetNumberOfIds();

  	// Current normal.
  	normals->GetTuple(i, n1);


  	for (j = 0; j < nearestPtCount; ++j){
  		combined->GetPoint(nearestPtIDs->GetId(j), b);

  		normals->GetTuple(nearestPtIDs->GetId(j), n2);

  		pt2ptDistSq = vtkMath::Distance2BetweenPoints(a, b);

  		val = constA * exp(constB * pt2ptDistSq);
  		val *= vtkMath::Dot(n1, n2);

  		totalDist += val;

  	}

  }

  totalDist = sqrt(totalDist);

  cout << totalDist << endl;

  return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
