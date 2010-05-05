/*=========================================================================

  Library   : packages
  Module    : $RCSfile: shreg.cc,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2000-2001
  Purpose   :
  Date      : $Date: 2004/12/17 16:25:37 $
  Version   : $Revision: 1.7 $
  Changes   : $Locker:  $
	      $Log: shreg.cc,v $
	      Revision 1.7  2004/12/17 16:25:37  dp1
	      Compilation with the new project
	
	      Revision 1.6  2004/12/16 16:58:47  dr
	      Introduced new transformation class hierachy
	
	      Revision 1.5  2004/11/15 11:40:03  tp500
	      Made program give more feedback on the data it is given.
	
	      Revision 1.3  2004/08/12 10:05:40  dr
	      Added compatibility for VTK 4.4 and higher

	      Revision 1.2  2004/07/22 10:54:12  dr
	      Added #ifdef HAS_VTK

	      Revision 1.1  2004/07/14 16:43:41  dr
	      Imported sources


=========================================================================*/

#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkFeatureEdges.h>
#include <vtkKDTreePointLocator.h>
#include <vtkGenericCell.h>
#include <vtkBMPReader.h>
#include <irtkLocator.h>
#include <vtkPolyDataWriter.h>

#include <irtkTransformation.h>

#include <irtkSurfaceRegistration.h>

//some global variables
char *_target_name = NULL, *_source_name = NULL;
char *dofin_name = NULL, *dofout_name = NULL;
char *fov_image_name = NULL;

void usage()
{
  cerr << "Usage: shreg [target] [source] [subdivisions] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-locator  value>         Locator: 0 = cell locator 1 = point locator 2 = kd-tree locator (default = point)" << endl;
  cerr << "<-dofin name>         Name of input file" << endl;
  cerr << "<-dofout name>        Name of output file" << endl;
  cerr << "<-ignoreedges>        Ignores edges in ICP (default OFF)" << endl;
  cerr << "<-epsilon>            Epsilon (default=0.01)" << endl;
  cerr << "<-numberofsteps>      Number of steps (default 5)" << endl;
  cerr << "<-stepsize>           Size of steps (default 5)" << endl;
  cerr << "<-iterations>         Number of iterations (default 100)" << endl;
  cerr << "<-subdivide>          Subdivide" << endl;
  cerr << "<-surfaces>           Number of surfaces in the data"<<endl;
  cerr << "<-ds>                 value" << endl;
  cerr << "<-fovImage>           Image to use when finding the fov bounds for dimensions of dofout." << endl;
  cerr << "<-fovDofin>           Use the details of dofin as the dimensions of dofout." << endl;
  exit(1);
}

void MarkBoundary(vtkPolyData *polydata)
{
  vtkFloatArray *scalars = vtkFloatArray::New();
  vtkFeatureEdges *edges = vtkFeatureEdges::New();
  edges->BoundaryEdgesOn();
  edges->FeatureEdgesOff();
  edges->ManifoldEdgesOff();
  edges->NonManifoldEdgesOff();
  edges->SetColoring(1);
  edges->SetInput(polydata);
  edges->Update();

  int size_before_filtering = polydata->GetNumberOfPoints();
  int size_after_filtering = edges->GetOutput()->GetNumberOfPoints();
  VTKFloat xyz[3];
  int closestPoint(0);

  vtkPointLocator *ploc = vtkPointLocator::New();
  ploc->SetDataSet(polydata);
  ploc->BuildLocator();

  scalars->SetNumberOfTuples(size_before_filtering);

  for (int counter = 0; counter < size_after_filtering; counter++){
    edges->GetOutput()->GetPoint(counter, xyz);
    closestPoint = ploc->FindClosestPoint(xyz);
    scalars->InsertTuple1(closestPoint, 0); //if edge set scalar for point to zero
  }

  for (int c = 0; c < size_before_filtering; c++){
    if (scalars->GetValue(c) != 0){
      scalars->InsertTuple1(c, 1);
    }
  }

  scalars->SetName("EDGEPOINTS");
  polydata->GetPointData()->SetScalars(scalars);
}

int main(int argc, char **argv)
{
  int i, elementsPerBucket, numberOfLevels;
  int locatorType, iterations, ok;
  float epsilon;
  double dx, dy, dz;
  bool ignoreEdges, subdivide;
  bool multipleSurfaces = false;
  int surfaces = 0;
  bool copyDofinDetails = false;

  irtkGreyImage fovImage;

  if (argc < 3){
    usage();
  }

  // Default parameters
  iterations = 100;
  elementsPerBucket = 5;
  locatorType = 1;
  epsilon = 0.01;
  ok = 0;
  ignoreEdges = False;
  subdivide = False;

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Parse filenames
  _target_name = argv[1];
  argv++;
  argc--;
  _source_name = argv[1];
  argv++;
  argc--;

  // Target pipeline
  cout << "Reading target ... " << _target_name;
  vtkPolyDataReader *target_reader = vtkPolyDataReader::New();
  target_reader->SetFileName(_target_name);
  target_reader->Modified();
  vtkPolyData *target = vtkPolyData::New();
  target = target_reader->GetOutput();
  target->Modified();
  target->Update();
  cout << "  (" << target->GetNumberOfPoints() << " points)" << endl;

  // Source pipeline
  cout << "Reading source ... " << _source_name;
  vtkPolyDataReader *source_reader = vtkPolyDataReader::New();
  source_reader->SetFileName(_source_name);
  source_reader->Modified();
  source_reader->Update();
  vtkPolyData *source = vtkPolyData::New();
  source = source_reader->GetOutput();
  source->Modified();
  cout << "  (" << source->GetNumberOfPoints() << " points)" << endl;

  // Number of subdivisions
  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)){
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-locator") == 0)){
      argc--;
      argv++;
      locatorType = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ignoreedges") == 0)){
      argc--;
      argv++;
      ignoreEdges = true;
      MarkBoundary(target);
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofout") == 0)){
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-epsilon") == 0)){
      argc--;
      argv++;
      epsilon = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-subdivide") == 0)){
      argc--;
      argv++;
      subdivide = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ds") == 0)){
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-surfaces") == 0)){
      argc--;
      argv++;
      multipleSurfaces = true;
      ok = True;      
    }
    if ((ok == False) && (strcmp(argv[1], "-fovImage") == 0)){
      argc--;
      argv++;
      fov_image_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-fovDofin") == 0)){
      argc--;
      argv++;
      copyDofinDetails = true;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkLocator *locator = new irtkLocator;
  locator->SelectLocatorType(locatorType);

  irtkSurfaceSimilarityNDDistanceMetric *metric = new irtkSurfaceSimilarityNDDistanceMetric;
  if (ignoreEdges) {
    metric->IgnoreEdges();
  }

  irtkSurfaceAffineRegistration *registration;

  if ( multipleSurfaces == true){
    registration  =    new irtkMSurfaceFreeFormRegistrationApproximation;
  }
  else
    registration =  new irtkSurfaceFreeFormRegistrationApproximation;
  

  registration->SetInput(target, source);
  registration->SetNumberOfIterations(iterations);
  registration->SetEpsilon(epsilon);
  registration->SetMetric(metric);
  registration->SetLocator(locator);

  // Create initial multi-level free-form deformation
  irtkMultiLevelBSplineFreeFormTransformation *mffd = NULL;

  // Read transformation
  if (dofin_name != NULL){
    if (irtkRigidTransformation::CheckHeader(dofin_name) == True){
      mffd = new irtkMultiLevelBSplineFreeFormTransformation;
      mffd->irtkRigidTransformation::Read(dofin_name);
    } else {
      if (irtkAffineTransformation::CheckHeader(dofin_name) == True){
	mffd = new irtkMultiLevelBSplineFreeFormTransformation;
	mffd->irtkAffineTransformation::Read(dofin_name);
      } else {
	if (irtkMultiLevelFreeFormTransformation::CheckHeader(dofin_name) == True){
	mffd = new irtkMultiLevelBSplineFreeFormTransformation(dofin_name);
	} else {
	  cerr << "Input transformation is not of type rigid, affine " << endl;
 	  cerr << "or multi-level free form deformation" << endl;
	  exit(1);
	}
      }
    }
  } else {
    // Otherwise use identity transformation to start
    mffd = new irtkMultiLevelBSplineFreeFormTransformation;
  }

  double start_box[3], end_box[3];
  VTKFloat *loc;
  double xaxis[3] = {1, 0, 0};
  double yaxis[3] = {0, 1, 0};

  if (fov_image_name != NULL){
    // Use the image to get a bounding box for the lattice.
    fovImage.Read(fov_image_name);

    start_box[0] = 0;
    start_box[1] = 0;
    start_box[2] = 0;
    end_box[0]   = fovImage.GetX()-1;
    end_box[1]   = fovImage.GetY()-1;
    end_box[2]   = fovImage.GetZ()-1;

    fovImage.ImageToWorld(start_box[0], start_box[1], start_box[2]);
    fovImage.ImageToWorld(end_box[0], end_box[1], end_box[2]);
    fovImage.GetOrientation(xaxis, yaxis);

  } else if (copyDofinDetails == true){
    // Generate a dof with the dimensions as the dofin.
    if (dofin_name == NULL){
      cout << "No dofin given from which to copy details." << endl;
      exit(1);
    }
    if (mffd->NumberOfLevels() > 1){
      cout << "Can only copy details from a single level FFD." << endl;
      exit(1);
    }
    irtkBSplineFreeFormTransformation *ffd = mffd->GetLocalTransformation(0);
    ffd->BoundingBox(start_box[0], start_box[1], start_box[2], end_box[0], end_box[1], end_box[2]);
    ffd->GetSpacing(dx, dy, dz);
    ffd->GetOrientation(xaxis, yaxis);
  } else {
    // Get bounding box of data
    loc = source->GetBounds();

    start_box[0] = *loc     - dx;
    start_box[1] = *(loc+2) - dy;
    start_box[2] = *(loc+4) - dz;
    end_box[0]   = *(loc+1) + dx;
    end_box[1]   = *(loc+3) + dy;
    end_box[2]   = *(loc+5) + dz;
  }

  cout << "Source Dimensions --  " << abs(start_box[0]-end_box[0]) << " x " << abs(start_box[1]-end_box[1]) <<
  					" x " << abs(start_box[2]-end_box[2]) << endl;
  cout << "Creating Grid     --  " << abs(round((start_box[0]-end_box[0])/dx))+1 << " x " << abs(round((start_box[1]-end_box[1])/dy))+1 <<
  					" x " << abs(round((start_box[2]-end_box[2])/dz))+1 << endl;
  cout << "Grid spacing      --  " << dx << "mm" << "  " << dy << "mm" << "  " << dz << "mm" << endl;
  // Create transformation
  irtkBSplineFreeFormTransformation *affd =  new
    irtkBSplineFreeFormTransformation(start_box[0], start_box[1], start_box[2],
				     end_box[0], end_box[1], end_box[2], dx, dy, dz, xaxis, yaxis);

  // Add ffd
  mffd->PushLocalTransformation(affd);

  for (i = 0; i < numberOfLevels-1; i++){

   // Set up registration and run
    registration->SetOutput(mffd);
    registration->Run();

    if (subdivide == False){
      // Add transformation
      dx = dx/2.0;
      dy = dy/2.0;
      dz = dz/2.0;
      irtkBSplineFreeFormTransformation *affd =  new
	irtkBSplineFreeFormTransformation(start_box[0], start_box[1], start_box[2],
					 end_box[0], end_box[1], end_box[2], dx, dy, dz);
      mffd->PushLocalTransformation(affd);
    } else {
     // Extract current transformation level
      irtkBSplineFreeFormTransformation *affd = 
	(irtkBSplineFreeFormTransformation *)mffd->GetLocalTransformation(0);

      // Subdivide
      affd->Subdivide();
    }
  }

  // Set up registration and run
  registration->SetOutput(mffd);
  registration->Run();

  // Write rigid transformation
  if (dofout_name != NULL){
    mffd->Write(dofout_name);
  }
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
