#include <irtkRegistrationMulti.h>

#define MAX_LEVELS 100

// Default filenames
char **source_names;
char **target_names;
char *dofin_name;
char *dofout_name[MAX_LEVELS];
char *parameter_name[MAX_LEVELS];

int paddingValue, numberOfLevels;

#ifdef HAS_VTK
extern Bool interactiveVTK;
extern Bool displayVTK;
extern Bool firstVTK;
#endif

void usage()
{
  cerr << "Usage: hregmulti [M] [N] [target1] ... [source1] ...\n"<< endl;
  cerr << "                 [subdivisions] <options> \n" << endl;

  cerr << "Register N source images to M target images."<< endl;
  cerr << "N.B. Currently only implemented for M = N (!)"<< endl;
  cerr << "Each target / source pair is considered a channel and "<< endl;
  cerr << "has its own similarity metric."<< endl;

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

void padding(irtkGreyImage *images, int nImages, 
             irtkMultiLevelFreeFormTransformation *mffd)
{
  int i, j, k, x, y, z, x1, y1, z1, x2, y2, z2, ok, index;
  int n;

  // Extract current transformation level
  irtkFreeFormTransformation *ffd = mffd->GetLocalTransformation(0);  

  // Calculate number of active and passive control points
  for (i = 0; i < ffd->GetX(); i++){
    for (j = 0; j < ffd->GetY(); j++){
      for (k = 0; k < ffd->GetZ(); k++){
	// Convert control points to index
	index = ffd->LatticeToIndex(i, j, k);

	// Calculate bounding box of control point in voxels
	ffd->BoundingBox(&images[0], index, x1, y1, z1, x2, y2, z2);

	ok = False;
	for (z = z1; z <= z2; z++){
	  for (y = y1; y <= y2; y++){
	    for (x = x1; x <= x2; x++){
              for (n = 0; n < nImages; ++n){
                if (images[n](x, y, z) > paddingValue){
		  ok = True;
	        }
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

void mreg(irtkGreyImage *targets, irtkGreyImage *sources, int nTargets, int nSources,
	  irtkMultiLevelBSplineFreeFormTransformation *mffd, int i)
{
  int count;
  irtkFreeFormRegistrationMulti registration;
  
  // Set input and output for the registration filter
  registration.SetNumberOfTargets(nTargets);
  registration.SetNumberOfSources(nSources);

  for (count = 0; count < nTargets; count++){
    registration.SetTarget(count, &targets[count]);
  }
  for (count = 0; count < nSources; count++){
    registration.SetSource(count, &sources[count]);
  }

  registration.SetOutput(mffd);

  // Set padding value
  if (paddingValue > MIN_GREY){
    registration.SetTargetPadding(paddingValue);
  }

  // Read default parameter
  registration.Read(parameter_name[i]);

  // Run registration filter
  registration.Run();
  
  // Write the final transformation estimate
  if (dofout_name[i] != NULL){
    mffd->Write(dofout_name[i]);
  }
  
  // Extract current transformation level
  irtkBSplineFreeFormTransformation *affd = mffd->GetLocalTransformation(0);  
  affd->Subdivide();
}

int main(int argc, char **argv)
{
  int i, ok, no_areg;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2; 
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2; 
  double x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis[3], yaxis[3];
  int nSources, nTargets;
  irtkAffineTransformation transformation;

  // Check command line
  if (argc < 6){
    usage();
  }

  // How many target and source images?
  nTargets = atoi(argv[1]);
  argc--;
  argv++;
  nSources  = atoi(argv[1]);
  argc--;
  argv++;

  // Currently only implemented for equal numbers of source and target
  // images.
  if (nSources != nTargets){
    usage();
  }

  // Parse source and target images
  target_names = new char*[nTargets];
  for (i = 0; i < nTargets; i++){
    target_names[i] = argv[1];
    argc--;
    argv++;
  }

  source_names = new char*[nSources];
  for (i = 0; i <nSources; i++){
    source_names[i] = argv[1];
    argc--;
    argv++;
  }

  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Read target images
  cout << "Reading targets ... "; cout.flush();
  irtkGreyImage *targets = new irtkGreyImage[nTargets];
  for (i = 0; i < nTargets; i++){
    targets[i].Read(target_names[i]);
  }
  cout << "done" << endl;
  // Read source images
  cout << "Reading source ... "; cout.flush();
  irtkGreyImage *sources = new irtkGreyImage[nSources];
  for (i = 0; i < nSources; i++){
    sources[i].Read(source_names[i]);
  }
  cout << "done" << endl;

  // Fix padding 
  paddingValue = MIN_GREY;

  // Fix ROI
  target_x1 = 0;
  target_y1 = 0;
  target_z1 = 0;
  target_x2 = targets[0].GetX();
  target_y2 = targets[0].GetY();
  target_z2 = targets[0].GetZ();
  source_x1 = 0;
  source_y1 = 0;
  source_z1 = 0;
  source_x2 = sources[0].GetX();
  source_y2 = sources[0].GetY();
  source_z2 = sources[0].GetZ();

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Parse remaining parameters
  no_areg = False;
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
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // If there is an region of interest, use it
  if ((target_x1 != 0) || (target_x2 != targets[0].GetX()) ||
      (target_y1 != 0) || (target_y2 != targets[0].GetY()) || 
      (target_z1 != 0) || (target_z2 != targets[0].GetZ())){
    for (i = 0; i < nTargets; i++){
      targets[i] = targets[i].GetRegion(target_x1, target_y1, target_z1, 
			                target_x2, target_y2, target_z2);
    }
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != sources[0].GetX()) ||
      (source_y1 != 0) || (source_y2 != sources[0].GetY()) || 
      (source_z1 != 0) || (source_z2 != sources[0].GetZ())){
    for (i = 0; i < nSources; i++){
      sources[i] = sources[i].GetRegion(source_x1, source_y1, source_z1, 
			                source_x2, source_y2, source_z2);
    }
  }

  // Print some information
  for (i = 0; i < numberOfLevels; i++){
    cout << "Performing non-rigid registration using " 
	 << parameter_name[i] << " writing result to " << dofout_name[i] 
	 << endl;
  }

  // Create transformation
  irtkMultiLevelBSplineFreeFormTransformation *mffd = 
    new irtkMultiLevelBSplineFreeFormTransformation;

  // Read affine transformation
  if (dofin_name != NULL){
    mffd->irtkAffineTransformation::Read(dofin_name);
  }

  // Create ffd
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = targets[0].GetX()-1;
  y2 = targets[0].GetY()-1;
  z2 = targets[0].GetZ()-1;
  targets[0].ImageToWorld(x1, y1, z1);
  targets[0].ImageToWorld(x2, y2, z2);
  targets[0].GetOrientation(xaxis, yaxis);

  irtkBSplineFreeFormTransformation *affd = new 
    irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis);

  // Add ffd
  mffd->PushLocalTransformation(affd);

  i = 0;
  while (i < numberOfLevels){
    // Perform padding
    padding(targets, nTargets, mffd);

    // Perform non-rigid registration
    mreg(targets, sources, nTargets, nSources, mffd, i);

    // Next level
    i++;
  }

  delete []target_names;
  delete []source_names;
  delete []targets;
  delete []sources;

}


