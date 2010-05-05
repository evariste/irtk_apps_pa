/*=========================================================================
 
  Library   : packages/applications
  Module    : $RCSfile: hreg.cc,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2000-2001
  Purpose   :
  Date      : $Date: 2006/11/16 10:17:49 $
  Version   : $Revision: 1.12 $
  Changes   : $Locker:  $
	      $Log: hreg.cc,v $
	      Revision 1.12  2006/11/16 10:17:49  dr
	      Added support for new parameter file format
	
	      Revision 1.11  2006/03/03 11:00:54  dr
	      Removed comments
	
	      Revision 1.10  2006/03/03 11:00:33  dr
	      Fixed bug so that subdivided FFD is pushed back on stack
	
	      Revision 1.9  2005/12/23 12:39:45  dr
	      Remove BSpline and Linear FFDs
	
	      Revision 1.8  2005/12/22 15:20:07  dr
	      Updated to reflect changes to transformation classes
	
	      Revision 1.7  2005/06/14 13:02:42  dr
	      Added MPI support
	
	      Revision 1.6  2005/05/11 17:09:08  dr
	      Changed name of access method from Put to Set
	
	      Revision 1.5  2005/05/04 12:29:51  dr
	      Add support to read non-rigid transformations via -dofin
	
	      Revision 1.4  2004/12/16 16:58:46  dr
	      Introduced new transformation class hierachy
	
	      Revision 1.3  2003/10/09 17:25:47  dr
	      Added text to help messages
	
	      Revision 1.2  2003/04/17 14:10:51  dr
	      Merged branch
	
	      Revision 1.1.2.6  2003/04/03 17:21:17  dr
	      Derived irtkMFreeFormTransformation from
	      irtkAffineTransformation
	
	      Revision 1.1.2.5  2003/01/03 17:16:30  dr
	      Added support for images with arbitrary orientations
	
	      Revision 1.1.2.4  2002/10/01 15:53:13  dr
	      Fixed bug when calculating bounding box of control points
	      for padding
	
	      Revision 1.1.2.3  2002/07/30 15:27:14  dr
	      Removed options for affine registration before non-rigid
	      registration. To perform affine registration use areg.cc
	      instead
	
	      Revision 1.1.2.2  2002/03/26 18:36:22  dr
	      Added support for different types of affine registration
	      (including no affine registration)
	
	      Revision 1.1.2.1  2001/11/16 18:13:09  dr
	      Imported sources
	

=========================================================================*/

#include <irtkRegistration.h>

#define MAX_LEVELS 100

#ifdef HAS_MPI

#include <mpi.h>

int mpi_nprocs = 1, mpi_rank;

#endif

// Default filenames
char *source_name;
char *target_name;
char *dofin_name;
char *dofout_name[MAX_LEVELS];
char *parameter_name[MAX_LEVELS];

int paddingValue, numberOfLevels, debug;

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
  cerr << "Usage: hreg [target] [source] [subdivisions] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parameter files>   Parameter files" << endl;
  cerr << "<-dofout    files>   Final transformation estimates" << endl;
  cerr << "<-dofin     file>    Initial transformation estimate" << endl;
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
  cerr << "Note that you need to specify the number of times the FFD " << endl;
  cerr << "should be subdivided. You need to specify a filename for " << endl;
  cerr << "the parameters and final transformation estimate for each" << endl;
  cerr << "subdivision." << endl;
  exit(1);
}

void padding(irtkGreyImage image, irtkMultiLevelFreeFormTransformation *mffd)
{
  int i, j, k, x, y, z, x1, y1, z1, x2, y2, z2, ok, index;

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

	ok = False;
	for (z = z1; z <= z2; z++){
	  for (y = y1; y <= y2; y++){
	    for (x = x1; x <= x2; x++){
	      if (image(x, y, z) > paddingValue){
		ok = True;
	      }
	    }
	  }
	}
	if (ok == False){
	  ffd->PutStatus(i, j, k, _Passive);
	}
      }
    }
  }
}

void mreg(irtkGreyImage target, irtkGreyImage source, 
	  irtkMultiLevelFreeFormTransformation *mffd, int i, int debug)
{
  irtkTempRegistration registration;
  
  // Set input and output for the registration filter
  registration.SetInput(&target, &source);
  registration.SetOutput(mffd);

  // Set padding value
  if (paddingValue > MIN_GREY){
    registration.SetTargetPadding(paddingValue);
  }

  // Debug flag
  registration.SetDebugFlag(debug);

  // Read default parameter
  registration.irtkImageRegistration::Read(parameter_name[i]);

  if (debug == True){
    registration.SetDebugFlag(True);
  }

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
  irtkFreeFormTransformation *affd = mffd->PopLocalTransformation();
  affd->Subdivide();
  mffd->PushLocalTransformation(affd);
}

int main(int argc, char **argv)
{
  int i, ok, no_areg;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2; 
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2; 
  double dx, dy, dz;
  irtkAffineTransformation transformation;
  int debug = False;

  // Check command line
  if (argc < 3){
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
  source_name = argv[1];
  argc--;
  argv++;
  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;
  // Read source image
  cout << "Reading source ... "; cout.flush();
  irtkGreyImage source(source_name);
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
  no_areg = False;
  debug   = False;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofout") == 0)){
      argc--;
      argv++;
      for (i = 0; i < numberOfLevels; i++){
	dofout_name[i] = argv[1];
	argc--;
	argv++;
      }
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)){
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)){
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)){
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)){
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)){
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)){
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx1") == 0)){
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx2") == 0)){
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty1") == 0)){
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty2") == 0)){
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz1") == 0)){
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz2") == 0)){
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx1") == 0)){
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx2") == 0)){
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy1") == 0)){
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy2") == 0)){
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz1") == 0)){
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz2") == 0)){
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      paddingValue = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dx") == 0)){
      argc--;
      argv++;
      dx = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dy") == 0)){
      argc--;
      argv++;
      dy = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dz") == 0)){
      argc--;
      argv++;
      dz = atof(argv[1]);
      argc--;
      argv++;
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
    if ((ok == False) && (strcmp(argv[1], "-parameter") == 0)){
      argc--;
      argv++;
      ok = True;
      for (i = 0; i < numberOfLevels; i++){
	parameter_name[i] = argv[1];
	argc--;
	argv++;
      }
    }
    if ((ok == False) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      debug = True;
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      ok = True;
      debug = True;
    }
    if (ok == False){
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

  // Read transformation
  if (dofin_name != NULL){
    irtkTransformation *transform = irtkTransformation::New(dofin_name);
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
    delete transform;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // Create ffd
  irtkBSplineFreeFormTransformation *affd = new irtkBSplineFreeFormTransformation(target, dx, dy, dz);

  // Add ffd
  mffd->PushLocalTransformation(affd);

  i = 0;
  while (i < numberOfLevels){
    // Perform padding
    padding(target, mffd);

    // Perform non-rigid registration
    mreg(target, source, mffd, i, debug);

    // Next level
    i++;
  }

#ifdef HAS_MPI
  cout << "Processor #" << mpi_rank << " spent " << (MPI_Wtime()-wtime) << "s" << endl;
  MPI_Finalize();
#endif

}


