#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#include <irtkRegistration.h>


#define MAX_LEVELS 100

#ifdef HAS_MPI

#include <mpi.h>

int mpi_nprocs = 1, mpi_rank;

#endif

// Default filenames
char *source_name;
char *source_surface_name;
char *target_name;
char *target_landmarks_name;
char *dofin_name;
char *dofout_name[MAX_LEVELS];
char *parameter_name[MAX_LEVELS];

int paddingValue, numberOfLevels;

#ifdef HAS_VTK
extern Bool interactiveVTK;
extern Bool displayVTK;
extern Bool firstVTK;
#endif

#ifdef HAS_MPI
  double wtime;
#endif

void usage()
{
  cerr << "Usage: constrainedhreg [target] [targetLandmarks] [source] [sourceSurface] [subdivisions] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parameter files>   Parameter files" << endl;
  cerr << "<-dofout    files>   Final transformation estimates" << endl;
  cerr << "<-dofin     file>    Initial transformation estimate" << endl;
  cerr << "<-ffdestimate>       Input transformation (dofin) must exist and should " << endl;
  cerr << "                     be of type mffd with a single level. The contained" << endl;
  cerr << "                     ffd will be used as the initial starting estimate" << endl;
  cerr << "                     for the control points." << endl;

  cerr << "<-Rx1 value>         Region of interest in both images" << endl;
  cerr << "<-Ry1 value>         Region of interest in both images" << endl;
  cerr << "<-Rz1 value>         Region of interest in both images" << endl;
  cerr << "<-Rx2 value>         Region of interest in both images" << endl;
  cerr << "<-Ry2 value>         Region of interest in both images" << endl;
  cerr << "<-Rz2 value>         Region of interest in both images" << endl;
  cerr << "<-Tx1 value>         Region of interest in target image" << endl;
  cerr << "<-Ty1 value>         Region of interest in target image" << endl;
  cerr << "<-Tz1 value>         Region of interest in target image" << endl;
  cerr << "<-Tx2 value>         Region of interest in target image" << endl;
  cerr << "<-Ty2 value>         Region of interest in target image" << endl;
  cerr << "<-Tz2 value>         Region of interest in target image" << endl;
  cerr << "<-Sx1 value>         Region of interest in source image" << endl;
  cerr << "<-Sy1 value>         Region of interest in source image" << endl;
  cerr << "<-Sz1 value>         Region of interest in source image" << endl;
  cerr << "<-Sx2 value>         Region of interest in source image" << endl;
  cerr << "<-Sy2 value>         Region of interest in source image" << endl;
  cerr << "<-Sz2 value>         Region of interest in source image" << endl;
  cerr << "<-Tp  value>         Padding value in target image" << endl;
  cerr << "<-ds  value>         Initial control point spacing" << endl;
  cerr << "<-dx  value>         Initial control point spacing for x" << endl;
  cerr << "<-dy  value>         Initial control point spacing for y" << endl;
  cerr << "<-dz  value>         Initial control point spacing for z" << endl;
  cerr << "<-debug>             Enable debugging information" << endl;

  cerr << "<-gradient>          Optimization by gradient descent" << endl;
  cerr << "<-conjugate>         Optimization by conjugate gradient descent" << endl;


  cerr << "Note that you need to specify the number of times the FFD " << endl;
  cerr << "should be subdivided. You need to specify a filename for " << endl;
  cerr << "the parameters and final transformation estimate for each" << endl;
  cerr << "subdivision." << endl;
  exit(1);
}

void padding(irtkGreyImage image, irtkMultiLevelFreeFormTransformation *mffd)
{
  int i, j, k, x, y, z, x1, y1, z1, x2, y2, z2, index;
  bool ok;

  // Extract current transformation level
  irtkFreeFormTransformation *ffd = mffd->GetLocalTransformation(0);  

  // Calculate number of active and passive control points
  for (i = 0; i < ffd->GetX(); i++){
    for (j = 0; j < ffd->GetY(); j++){
      for (k = 0; k < ffd->GetZ(); k++){
	// Convert control points to index
	index = ffd->LatticeToIndex(i, j, k);

	// Calculate bounding box of control point in voxels
	ffd->BoundingBox(&image, index, x1, y1, z1, x2, y2, z2);

	ok = false;
	for (z = z1; z <= z2; z++){
	  for (y = y1; y <= y2; y++){
	    for (x = x1; x <= x2; x++){
	      if (image(x, y, z) > paddingValue){
		ok = true;
	      }
	    }
	  }
	}
	if (ok == false){
	  ffd->PutStatus(i, j, k, _Passive);
	}
      }
    }
  }
}

void mreg(irtkGreyImage target, vtkPolyData *targetLandmarks, irtkGreyImage source, 
	  vtkPolyData *sourceSurface, irtkMultiLevelFreeFormTransformation *mffd, int i, irtkOptimizationMethod optimizationMethod)
{
  irtkTempRegistration registration;
  
  // Set input and output for the registration filter
  registration.SetInput(&target, &source);
  registration.SetOutput(mffd);

  registration.SetTargetLandmarks(targetLandmarks);
  registration.SetSourceSurface(sourceSurface);

  registration.SetOptimizationMethod(optimizationMethod);

  // Set padding value
  if (paddingValue > MIN_GREY){
    registration.SetTargetPadding(paddingValue);
  }

  // Read default parameter
  registration.Read(parameter_name[i]);

  // Run registration filter
  registration.Run();
  
  // Write the final transformation estimate
#ifdef HAS_MPI
  if (mpi_rank == 0)
#endif
  {
    if (dofout_name[i] != NULL){
      mffd->irtkTransformation::Write(dofout_name[i]);
    }
  }
  
  // Extract current transformation level
  irtkFreeFormTransformation *affd = mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);  
  affd->Subdivide();
}

int main(int argc, char **argv)
{
  int i, ok, no_areg;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2; 
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2; 
  double x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis[3], yaxis[3];

  irtkOptimizationMethod optimizationMethod = GradientDescent;

  irtkAffineTransformation transformation;

  int use_FFD_as_starting_estimate = false;

  // Check command line
  if (argc < 6){
    usage();
  }

#ifdef HAS_MPI
  MPI_Init(&argc, &argv);
  wtime = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  target_landmarks_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  source_surface_name = argv[1];
  argc--;
  argv++;
  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;

  // Read target landmarks
  cout << "Reading target landmarks ... "; cout.flush();
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(target_landmarks_name);
  reader->Update();
  vtkPolyData* targetLandmarks = vtkPolyData::New();
  targetLandmarks = reader->GetOutput();
  cout << "done" << endl;

  // Read source image
  cout << "Reading source ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

  // Read source surface
  cout << "Reading source surface ... "; cout.flush();
  vtkPolyDataReader *reader2 = vtkPolyDataReader::New();
  reader2->SetFileName(source_surface_name);
  reader2->Update();
  vtkPolyData* sourceSurface = vtkPolyData::New();
  sourceSurface = reader2->GetOutput();
  cout << "done" << endl;

  // Fix padding 
  paddingValue = MIN_GREY;

  // Fix ROI
  target_x1 = 0;
  target_y1 = 0;
  target_z1 = 0;
  target_x2 = target.GetX();
  target_y2 = target.GetY();
  target_z2 = target.GetZ();
  source_x1 = 0;
  source_y1 = 0;
  source_z1 = 0;
  source_x2 = source.GetX();
  source_y2 = source.GetY();
  source_z2 = source.GetZ();

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Parse remaining parameters
  no_areg = false;
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ffdestimate") == 0)){
      argc--;
      argv++;
      use_FFD_as_starting_estimate = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)){
      argc--;
      argv++;
      for (i = 0; i < numberOfLevels; i++){
	dofout_name[i] = argv[1];
	argc--;
	argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)){
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)){
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)){
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)){
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)){
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)){
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)){
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)){
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)){
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)){
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)){
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)){
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sx1") == 0)){
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sx2") == 0)){
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sy1") == 0)){
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sy2") == 0)){
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sz1") == 0)){
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sz2") == 0)){
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      paddingValue = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dx") == 0)){
      argc--;
      argv++;
      dx = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dy") == 0)){
      argc--;
      argv++;
      dy = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dz") == 0)){
      argc--;
      argv++;
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)){
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-parameter") == 0)){
      argc--;
      argv++;
      ok = true;
      for (i = 0; i < numberOfLevels; i++){
	parameter_name[i] = argv[1];
	argc--;
	argv++;
      }
    }
    if ((ok == false) && (strcmp(argv[1], "-conjugate") == 0)){
      argc--;
      argv++;
      ok = true;
      optimizationMethod = ConjugateGradientDescentWithConstraints;
    }
    if ((ok == false) && (strcmp(argv[1], "-gradient") == 0)){
      argc--;
      argv++;
      ok = true;
      optimizationMethod = GradientDescentWithConstraints;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // If there is an region of interest, use it
  if ((target_x1 != 0) || (target_x2 != target.GetX()) ||
      (target_y1 != 0) || (target_y2 != target.GetY()) || 
      (target_z1 != 0) || (target_z2 != target.GetZ())){
    target = target.GetRegion(target_x1, target_y1, target_z1, 
			      target_x2, target_y2, target_z2);
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != source.GetX()) ||
      (source_y1 != 0) || (source_y2 != source.GetY()) || 
      (source_z1 != 0) || (source_z2 != source.GetZ())){
    source = source.GetRegion(source_x1, source_y1, source_z1, 
			      source_x2, source_y2, source_z2);
  }

  // Print some information
  for (i = 0; i < numberOfLevels; i++){
    cout << "Performing non-rigid registration using " 
	 << parameter_name[i] << " writing result to " << dofout_name[i] 
	 << endl;
  }

  // Create transformation
  irtkMultiLevelFreeFormTransformation *mffd = 
    new irtkMultiLevelFreeFormTransformation;

  irtkTransformation *transform = NULL;

  // Read transformation
  if (dofin_name != NULL){
    transform = irtkTransformation::New(dofin_name);
    if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0){
      mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
    } else {
      if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0){
	mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
      } else {
	if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0){
	  mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
	} else {
	  cerr << "Input transformation is not of type rigid, affine " << endl;
 	  cerr << "or multi-level free form deformation" << endl;
	  exit(1);
	}
      }
    }
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  if (use_FFD_as_starting_estimate = true){
    // Check a dofin with exactly one level was given.
    if (mffd->NumberOfLevels() != 1){
      cerr << "-ffdestimate flag given but input transformation " << endl;
      cerr << "does not have exactly one level." << endl;
      exit(1);
    }
    cerr << "Using FFD to initialise control point values." << endl;
  } else {
    // Create ffd
    x1 = 0;
    y1 = 0;
    z1 = 0;
    x2 = target.GetX()-1;
    y2 = target.GetY()-1;
    z2 = target.GetZ()-1;
    target.ImageToWorld(x1, y1, z1);
    target.ImageToWorld(x2, y2, z2);
    target.GetOrientation(xaxis, yaxis);

    irtkBSplineFreeFormTransformation *affd = new 
      irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis);

    // Add ffd
    mffd->PushLocalTransformation(affd);
  } 

  i = 0;
  while (i < numberOfLevels){
    // Perform padding
    padding(target, mffd);

    // Perform non-rigid registration
    mreg(target, targetLandmarks, source, sourceSurface, mffd, i, optimizationMethod);

    // Next level
    i++;
  }

  if (transform != NULL)
    delete transform;


#ifdef HAS_MPI
  cout << "Processor #" << mpi_rank << " spent " << (MPI_Wtime()-wtime) << "s" << endl;
  MPI_Finalize();
#endif

}



#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
