#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

// Default filenames
char *source_name = NULL, *output_name = NULL, *target_name = NULL;

char **dof_name  = NULL;

#define MAX_DOFS 50

void usage()
{
  cerr << " Usage: transformationN [source] [output] <options>\n" << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " " << endl;
  cerr << " <-target image>    ." << endl;
  cerr << " " << endl;
  cerr << " <-dofin file>      Transformation file. Multiple transformations can be given by" << endl;
  cerr << "                    repeatedly using this flag. They are processed in order they passed." << endl;
  cerr << " <-dofin_i file>    A transformation whose inverse should be applied." << endl;
  cerr << " " << endl;
  cerr << " <-Tp  value>       Padding value in target." << endl;
  cerr << " <-Sp  value>       Padding value in source (default 0)" << endl;
  cerr << " " << endl;
  cerr << " <-linear>          Linear interpolation" << endl;
  cerr << " <-bspline>         B-spline interpolation" << endl;
  cerr << " " << endl;
  cerr << " E.g." << endl;
  cerr << " " << endl;
  cerr << "   transformationN src.nii.gz out.nii.gz -target tgt.nii.gz -dofin a.dof -dofin_i b.dof -dofin c.dof" << endl;
  cerr << " " << endl;
  cerr << " reslices the source using the composed transformation c b^-1 a(x) applied to each target location x." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation **transformation = NULL;
  irtkInterpolateImageFunction *interpolator = NULL;

  bool ok;
  int i, j, k;
  double x, y, z;

  double x1, y1, z1, x2, y2, z2, target_padding, source_padding, val;
  int m, noOfDofs;

  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse source and target images
  source_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read source image
  irtkRealImage source(source_name);

  target_padding = -1.0 * FLT_MAX;
  source_padding = 0;

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

    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      target_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sp") == 0)) {
      argc--;
      argv++;
      source_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }


    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)){
      argc--;
      argv++;

      if (threeD){
      	interpolator = new irtkLinearInterpolateImageFunction;
      } else {
      	interpolator = new irtkLinearInterpolateImageFunction2D;
      	cerr << "Using 2D interpolator" << endl;
      }

      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)){
      argc--;
      argv++;

      if (threeD){
      	interpolator = new irtkBSplineInterpolateImageFunction;
      } else {
      	interpolator = new irtkBSplineInterpolateImageFunction2D;
      }

      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      target_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkRealImage target;

  // Set up a target,
  if (target_name != NULL){
    target.Read(target_name);
  } else {
    target.Read(source_name);
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



  // Fill histogram
  for (k = 0; k < target.GetZ(); k++){
    for (j = 0; j < target.GetY(); j++){
      for (i = 0; i < target.GetX(); i++){

      	irtkPoint p(i, j, k);

      	// Convert point into world coordinates
      	target.ImageToWorld(p);

      	// Transform point
      	for (m = 0; m < noOfDofs; m++){
      		if (invert[m]){
      			x = p._x;
      			y = p._y;
      			z = p._z;
        		transformation[m]->Inverse(x, y, z);
        		p._x = x;
        		p._y = y;
        		p._z = z;
      		} else {
      			transformation[m]->Transform(p);
      		}
      	}

      	// Convert point into image coordinates
      	source.WorldToImage(p);

      	// Put in the fall back value in case we've ended up outside the source
      	// image's field of view.
      	target.Put(i, j, k, source_padding);

      	// A bad thing might happen for the 2D case.
      	if (!threeD &&
      			(p._z > 0.5 || p._z < -0.5)){
      		cerr << "Transformed point outside plane of 2D source image." << endl;
      		exit(1);
      	}

      	// Conditions for leaving the fall back value.
      	if (threeD &&
      			(p._x <= x1 || p._x >= x2 ||
      			 p._y <= y1 || p._y >= y2 ||
      			 p._z <= z1 || p._z >= z2))
      		continue;

      	if (!threeD &&
      			(p._x <= x1 || p._x >= x2 ||
      			 p._y <= y1 || p._y >= y2))
      		continue;

      	// If we've got this far, we can safely interpolate in the
      	// source image.
      	val = interpolator->EvaluateInside(p._x, p._y, p._z);

      	target.Put(i, j, k, val);

      }
    }
  }

  target.Write(output_name);

  delete [] dof_name;
  delete [] invert;

}

