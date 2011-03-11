#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

// Default filenames
char *source_name = NULL, *dofOut_name = NULL;

char **dof_name  = NULL;

#define MAX_DOFS 50

void usage()
{
  cerr << " Usage: ffdcomposeN [source] [dofOut] <options>\n" << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " " << endl;
  cerr << " <-dofin file>      Transformation file. Multiple transformations can be given by" << endl;
  cerr << "                    repeatedly using this flag. They are processed in order they passed." << endl;
  cerr << " <-dofin_i file>    A transformation whose inverse should be applied." << endl;
  cerr << " " << endl;
  cerr << " <-linear>          Linear interpolation" << endl;
  cerr << " <-bspline>         B-spline interpolation" << endl;
  cerr << " " << endl;
  cerr << " E.g." << endl;
  cerr << " " << endl;
  cerr << "   ffdcomposeN out.dof.gz -dofin a.dof -dofin_i b.dof -dofin c.dof" << endl;
  cerr << " " << endl;
  cerr << " returns the composed transformation c b^-1 a(x) applied to each target location x." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation **transformation = NULL;
  irtkInterpolateImageFunction *interpolator = NULL;

  int i, j, k, m, xdim, ydim, zdim, noOfDofs, numberOfCPs, count;
  double x, y, z;
  double x1, y1, z1, x2, y2, z2;

  bool ok;

  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse dofOut_name
  source_name = argv[1];
  argc--;
  argv++;
  dofOut_name = argv[1];
  argc--;
  argv++;

  // Read source image
  irtkRealImage source(source_name);

  bool threeD = true;

  if (source.GetZ() == 1){
  	threeD = false;
  }

  // Fix number of dofs
  noOfDofs = 0;

  dof_name = new char*[MAX_DOFS];
  bool *invert = new bool[MAX_DOFS];

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      invert[noOfDofs]   = false;
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin_i") == 0)){
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      invert[noOfDofs]   = true;
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)){
      argc--;
      argv++;
      ok = true;

      if (threeD){
      	interpolator = new irtkLinearInterpolateImageFunction;
      } else {
      	interpolator = new irtkLinearInterpolateImageFunction2D;
      	cerr << "Using 2D interpolator" << endl;
      }
    }

    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)){
      argc--;
      argv++;
      ok = true;
 
      if (threeD){
      	interpolator = new irtkBSplineInterpolateImageFunction;
      } else {
      	interpolator = new irtkBSplineInterpolateImageFunction2D;
      }
   }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Create default interpolator if necessary
  if (interpolator == NULL){
  	if (threeD)
  		interpolator = new irtkNearestNeighborInterpolateImageFunction;
  	else
  		interpolator = new irtkNearestNeighborInterpolateImageFunction2D;
  }

  interpolator->SetInput(&source);
  interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  interpolator->Inside(x1, y1, z1, x2, y2, z2);

  if (noOfDofs == 0){
    noOfDofs = 1;
    transformation = new irtkTransformation*[noOfDofs];
    transformation[0] = new irtkRigidTransformation;
  } else {
    transformation = new irtkTransformation*[noOfDofs];
    for (m = 0; m < noOfDofs; m++){
      transformation[m] = irtkTransformation::New(dof_name[m]);
    }
  }

  irtkMultiLevelFreeFormTransformation *mffd_out = new irtkMultiLevelFreeFormTransformation;
  mffd_out->irtkTransformation::Read(dof_name[noOfDofs-1]);  // takes the last DOF, for now

  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *affd_out = dynamic_cast<irtkFreeFormTransformation3D *>(mffd_out->GetLocalTransformation(0));
  xdim = affd_out->GetX();
  ydim = affd_out->GetY();
  zdim = affd_out->GetZ();
  numberOfCPs = xdim * ydim * zdim;

  // Space to store the control point displacements.
  double *xdata = new double[numberOfCPs];
  double *ydata = new double[numberOfCPs];
  double *zdata = new double[numberOfCPs];

  count = 0;

  // Loop for each control point in the target
  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        x = i; // why? i not used later
        y = j; // why? j not used later
        z = k; // why? k not used later

        // Transform point from lattice coordinates to target coordinates
        affd_out->LatticeToWorld(x, y, z);

      	// Transform point
      	for (m = 0; m < noOfDofs; m++){
	  if (invert[m])
	    transformation[m]->Inverse(x, y, z);
	  else
	    transformation[m]->Transform(x,y,z);
      	}

        // New location of the control point.
        xdata[count] = x;
        ydata[count] = y;
        zdata[count] = z;

        ++count; 

      	// A bad thing might happen for the 2D case.
      	if (!threeD &&
      			(z > 0.5 || z < -0.5)){
      		cerr << "Transformed point outside plane of 2D source image." << endl;
      		exit(1);
      	}

      	// Conditions for leaving the fall back value.
      	if (threeD &&
      			(x <= x1 || x >= x2 ||
      			 y <= y1 || y >= y2 ||
      			 z <= z1 || z >= z2))
      		continue;

      	if (!threeD &&
      			(x <= x1 || x >= x2 ||
      			 y <= y1 || y >= y2))
      		continue;

      }
    }
  }

  // Interpolate the ffd and write dof
  affd_out->Interpolate(xdata, ydata, zdata);
  mffd_out->irtkTransformation::Write(dofOut_name);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;
  delete [] dof_name;
  delete [] invert;

}

