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
char *maskName = NULL;


void usage()
{
  cerr << " polydatasurfacedistance [inputA] [inputB]" << endl;
  cerr << " " << endl;
  cerr << " Estimate the 'currents' type distance between a pair" << endl;
  cerr << " of surfaces. See the Glaunes IPMI 2005 paper and " << endl;
  cerr << " Linh Ha, MICCAI, 2010." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " " << endl;
  cerr << " -mask name  : Name of scalars mask. Restrict distance summation to faces with positive values of the mask." << endl;
  cerr << "               Scalars with the same name must be present in both surfaces" << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  cerr << " " << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, jj;
  bool ok, useMask;
  double searchRadius;
  double pt1[3], pt2[3], pt3[3];
  double e1[3], e2[3], e3[3];
  double n1[3], n2[3];
  double volA, volB, rA, rB, sigmaKer, kernelFraction;
  int nearestPtCount;
  int noOfFacesA, noOfFacesB;
  double boundsA[6], boundsB[6];
  double gaussConst1, gaussConst2;
  double val;

  bool verbose = false;

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
    ok = false;
    if ((!ok) && (strcmp(argv[1], "-kernelFraction") == 0)) {
      argc--;
      argv++;
      kernelFraction = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-mask") == 0)) {
      argc--;
      argv++;
      maskName = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-verbose") == 0)) {
          argc--;
          argv++;
          verbose = true;
          ok = true;
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

  vtkFloatArray *maskA = vtkFloatArray::New();
  vtkFloatArray *maskB = vtkFloatArray::New();

  useMask = false;
  int ind;

  if (maskName != NULL){
    maskA = (vtkFloatArray*) triAout->GetPointData()->GetArray(maskName, ind);
    if (ind == -1 || maskA == NULL){
      cerr << "Scalars unavailable with name " << maskName << " in surface " << input_nameA << endl;
      exit(1);
    }
    maskB = (vtkFloatArray*) triBout->GetPointData()->GetArray(maskName, ind);
    if (ind == -1 || maskB == NULL){
      cerr << "Scalars unavailable with name " << maskName << " in surface " << input_nameB << endl;
      exit(1);
    }

    useMask =true;
  }

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


  vtkPoints *centresA = vtkPoints::New();
  centresA->SetNumberOfPoints(noOfFacesA);

  vtkPoints *centresB = vtkPoints::New();
  centresB->SetNumberOfPoints(noOfFacesB);


  vtkFloatArray *normalsA = vtkFloatArray::New();
  normalsA->SetNumberOfComponents(3);
  normalsA->SetNumberOfTuples(noOfFacesA);

  vtkFloatArray *normalsB = vtkFloatArray::New();
  normalsB->SetNumberOfComponents(3);
  normalsB->SetNumberOfTuples(noOfFacesB);


  vtkFloatArray *maskOutA = vtkFloatArray::New();
  maskOutA->SetNumberOfComponents(1);
  maskOutA->SetNumberOfTuples(noOfFacesA);

  vtkFloatArray *maskOutB = vtkFloatArray::New();
  maskOutB->SetNumberOfComponents(1);
  maskOutB->SetNumberOfTuples(noOfFacesB);




  vtkIdType nptsForFace = 0;
  vtkIdType *ptIdsForFace;
  vtkIdType ptIDa, ptIDb, ptIDc;


  facesA->InitTraversal();
  for (i = 0; i < noOfFacesA; ++i){

  	facesA->GetNextCell(nptsForFace, ptIdsForFace);

  	if (nptsForFace != 3){
  		cerr << "Error: faces must have three points each." << endl;
  		exit(0);
  	}

    ptIDa = ptIdsForFace[0];
    ptIDb = ptIdsForFace[1];
    ptIDc = ptIdsForFace[2];

  	triAout->GetPoint(ptIDa, pt1);
  	triAout->GetPoint(ptIDb, pt2);
  	triAout->GetPoint(ptIDc, pt3);

  	for (j = 0; j < 3; ++j){
  		e1[j] = pt2[j] - pt1[j];
  		e2[j] = pt3[j] - pt1[j];
  	}

  	vtkMath::Cross(e1, e2, e3);
  	normalsA->SetTuple3(i, e3[0], e3[1], e3[2]);


  	// Get centroid and set masking.
  	for (j = 0; j < 3; ++j){
  		pt1[j] += (pt2[j] + pt3[j]);
  		pt1[j] /= 3.0;
  	}

  	centresA->SetPoint(i, pt1);

    val = 0;
    if (useMask){
      if( maskA->GetTuple1(ptIDa) > 0 )
        val++;
      if( maskA->GetTuple1(ptIDb) > 0 )
        val++;
      if( maskA->GetTuple1(ptIDc) > 0 )
        val++;
    }

  	maskOutA->SetTuple1(i, (val > 1) ? 1 : 0);

  }


  facesB->InitTraversal();

  // Note where this loop starts and ends.
  for (i = 0; i < noOfFacesB; ++i){

  	facesB->GetNextCell(nptsForFace, ptIdsForFace);

  	if (nptsForFace != 3){
  	  		cerr << "Error: faces must have three points each." << endl;
  	  		exit(0);
  	}

    ptIDa = ptIdsForFace[0];
    ptIDb = ptIdsForFace[1];
    ptIDc = ptIdsForFace[2];

    triBout->GetPoint(ptIDa, pt1);
    triBout->GetPoint(ptIDb, pt2);
    triBout->GetPoint(ptIDc, pt3);

  	for(j = 0; j < 3; ++j){
  		e1[j] = pt2[j] - pt1[j];
  		e2[j] = pt3[j] - pt1[j];
  	}

  	// Cross product is skew symmetric for the second surface, we do e2 x e1 while we did
  	// e1 x e2 for the first surface - i.e.  the sense is negated.
  	vtkMath::Cross(e2, e1, e3);
  	normalsB->SetTuple3(i, e3[0], e3[1], e3[2]);

    // Get centroid and set masking.
  	for (j = 0; j < 3; ++j){
  		pt1[j] += (pt2[j] + pt3[j]);
  		pt1[j] /= 3.0;
  	}

  	centresB->SetPoint(i, pt1);

    val = 0;
    if (useMask){
      if( maskB->GetTuple1(ptIDa) > 0 )
        val++;
      if( maskB->GetTuple1(ptIDb) > 0 )
        val++;
      if( maskB->GetTuple1(ptIDc) > 0 )
        val++;
    }
    maskOutB->SetTuple1(i, (val > 1) ? 1 : 0);
  }

  normalsA->SetName("faceNormals");
  maskOutA->SetName("maskOut");

  normalsB->SetName("faceNormals");
  maskOutB->SetName("maskOut");


  vtkPolyData *currentA = vtkPolyData::New();
  currentA->SetPoints(centresA);
  currentA->GetPointData()->AddArray(normalsA);
  currentA->GetPointData()->AddArray(maskOutA);
  currentA->Update();

  vtkPolyData *currentB = vtkPolyData::New();
  currentB->SetPoints(centresB);
  currentB->GetPointData()->AddArray(normalsB);
  currentB->GetPointData()->AddArray(maskOutB);
  currentB->Update();


  vtkPointLocator *point_locatorA = vtkPointLocator::New();
  point_locatorA->SetNumberOfPointsPerBucket(5);
  point_locatorA->SetDataSet(currentA);
  point_locatorA->BuildLocator();

  vtkPointLocator *point_locatorB = vtkPointLocator::New();
  point_locatorB->SetNumberOfPointsPerBucket(5);
  point_locatorB->SetDataSet(currentB);
  point_locatorB->BuildLocator();



  // Need to establish a suitable radius for the kernel on centre
  // to centre distances.

  volA = (boundsA[1] - boundsA[0]) *
				 (boundsA[3] - boundsA[2]) *
  		   (boundsA[5] - boundsA[4]);

  volB = (boundsB[1] - boundsB[0]) *
  		   (boundsB[3] - boundsB[2]) *
  		   (boundsB[5] - boundsB[4]);

  // Assume points distributed on spheres.
  rA = pow(volA * 3.0 / 4.0 / M_PI, (1.0/3.0));
  rB = pow(volB * 3.0 / 4.0 / M_PI, (1.0/3.0));


  sigmaKer = 0.5 * (rA + rB) * kernelFraction;

  if (verbose){
    cout << "Estimated radii  : " << rA << " and " << rB << endl;
    cout << "Sigma for kernel : " << sigmaKer << endl;
  }

  // Constants for kernel
  gaussConst1 = 1 / sigmaKer / sqrt(2 * M_PI);
  gaussConst2 = -1.0 / sigmaKer / sigmaKer;

  searchRadius = 2.5 * sigmaKer;


  vtkIdList *nearestPtIDs = vtkIdList::New();

  double pt2ptDistSq, totalDist;
  double minusADotB, minusBDotA, aDotA, bDotB;

  // The distance we seek to measure and return.
  totalDist = 0.0;


  // inter surface distance:
  for (i = 0; i < noOfFacesA; ++i){

    if (useMask && maskOutA->GetTuple1(i) <= 0){
    	continue;
    }

  	currentA->GetPoint(i, pt1);

  	point_locatorB->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);

  	nearestPtCount = nearestPtIDs->GetNumberOfIds();

  	// Current normal.
  	normalsA->GetTuple(i, n1);

  	for (j = 0; j < nearestPtCount; ++j){
  		jj = nearestPtIDs->GetId(j);

  		currentB->GetPoint(jj, pt2);
  		normalsB->GetTuple(jj, n2);

  		pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);

  		val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
  		val *= vtkMath::Dot(n1, n2);

  		totalDist += val;

  	}
  }

  minusADotB = totalDist;
  totalDist = 0.0;

  for (i = 0; i < noOfFacesB; ++i){

    if (useMask && maskOutB->GetTuple1(i) <= 0){
    	continue;
    }

  	currentB->GetPoint(i, pt1);

  	point_locatorA->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);

  	nearestPtCount = nearestPtIDs->GetNumberOfIds();

  	// Current normal.
  	normalsB->GetTuple(i, n1);

  	for (j = 0; j < nearestPtCount; ++j){
  		jj = nearestPtIDs->GetId(j);

  		currentA->GetPoint(jj, pt2);
  		normalsA->GetTuple(jj, n2);

  		pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);

  		val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
  		val *= vtkMath::Dot(n1, n2);

  		totalDist += val;

  	}
  }

  minusBDotA = totalDist;
  totalDist = 0.0;

// Within dists for checking : to drop

//  for (i = 0; i < noOfFacesA; ++i){
//
//    if (useMask && maskOutA->GetTuple1(i) <= 0){
//    	continue;
//    }
//
//  	currentA->GetPoint(i, pt1);
//
//  	point_locatorA->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);
//
//  	nearestPtCount = nearestPtIDs->GetNumberOfIds();
//
//  	// Current normal.
//  	normalsA->GetTuple(i, n1);
//
//  	for (j = 0; j < nearestPtCount; ++j){
//  		jj = nearestPtIDs->GetId(j);
//
//  		currentA->GetPoint(jj, pt2);
//  		normalsA->GetTuple(jj, n2);
//
//  		pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);
//
//  		val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
//  		val *= vtkMath::Dot(n1, n2);
//
//  		totalDist += val;
//
//  	}
//  }
//
//  aDotA = totalDist;
//  totalDist = 0.0;
//
//  for (i = 0; i < noOfFacesB; ++i){
//
//    if (useMask && maskOutB->GetTuple1(i) <= 0){
//    	continue;
//    }
//
//  	currentB->GetPoint(i, pt1);
//
//  	point_locatorB->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);
//
//  	nearestPtCount = nearestPtIDs->GetNumberOfIds();
//
//  	// Current normal.
//  	normalsB->GetTuple(i, n1);
//
//  	for (j = 0; j < nearestPtCount; ++j){
//  		jj = nearestPtIDs->GetId(j);
//
//  		currentB->GetPoint(jj, pt2);
//  		normalsB->GetTuple(jj, n2);
//
//  		pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);
//
//  		val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
//  		val *= vtkMath::Dot(n1, n2);
//
//  		totalDist += val;
//
//  	}
//  }
//
//  bDotB = totalDist;


  cout << minusADotB << " " << minusBDotA << " " << sqrt(fabs(minusADotB + minusBDotA)) << endl;


//  ////////////////////////////////////////////////////////////////
//  tempScalars->SetName("tempScalars");
//  combined->GetPointData()->AddArray(tempScalars);
//  combined->Update();
//
//  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//  writer->SetInput(combined);
//  writer->SetFileTypeToBinary();
//  writer->SetFileName("temp.vtk");
//  writer->Write();
//  exit(0);
//  ////////////////////////////////////////////////////////////////


  return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
