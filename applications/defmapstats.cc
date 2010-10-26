
#ifdef HAS_VTK

#include <irtkRegistration.h>
#include <nr.h>

// Default filenames
char *target_name = NULL, *output_name = NULL, *mask_name = NULL;
char **dof_name  = NULL;

#define MAX_DOFS 10

#define MAX_PTS_PAREG 10000

void usage()
{
  cerr << " Usage: defmapstats [target] <options>" << endl;
  cerr << " " << endl;
  cerr << " Stats of distances that voxels in the target move under the effect" << endl;
  cerr << " of one or more transformations." << endl;
  cerr << " " << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " -dofin file      Transformation. Can be repeated to give multiple dofs" << endl;
  cerr << "                  (i.e. \"-dofin file1 -dofin file2 ...\") which are applied" << endl;
  cerr << "                  in order given (Max 10)." << endl;
  cerr << " -Tp  value       Padding value in target." << endl;
  cerr << " -mask file       Image mask to define a region of interest." << endl;
  cerr << " -removeGlobal    Estimate a global affine transformation based on a sampled" << endl;
  cerr << "                  subset of voxel locations and their transformed coordinates." << endl;
  cerr << "                  Remove the effect of this global transformation before " << endl;
  cerr << "                  estimating the distances." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation **transformation = NULL;
  bool ok;
  int i, j, k, regressAffine, squaredDistance;
  double Tp, val;
  int m, noOfDofs;
  int noOfPoints, ptID;
  double sumD, sumD2;
  double x, y, z, dx, dy, dz;

  int xdim, ydim, zdim;
  float *vals, *vals2, *L, *L2, sumL, sumL2;
  double alpha = 1;


  // Check command line
  if (argc < 2){
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;

  // Read target image
  irtkRealImage target(target_name);

  Tp = -1.0 * FLT_MAX;

  // Fix number of dofs
  noOfDofs = 0;

  dof_name = new char*[MAX_DOFS];

  regressAffine = false;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      Tp = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-removeGlobal") == 0)){
      argc--;
      argv++;
      regressAffine = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask_name = argv[1];
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
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (noOfDofs == 0){
    cerr << "No transformations specified: Using a single identity transformation!" << endl;
    noOfDofs = 1;
    transformation = new irtkTransformation*[noOfDofs];
    transformation[0] = new irtkRigidTransformation;
  } else {
    // Read dof(s)
    transformation = new irtkTransformation*[noOfDofs];
    for (m = 0; m < noOfDofs; m++){
      transformation[m] = irtkTransformation::New(dof_name[m]);
    }
  }

  // Setting up mask
  cout.flush();

  irtkGreyImage mask;

  if (mask_name == NULL){
    mask.Read(target_name);
    irtkGreyPixel *ptr2mask = mask.GetPointerToVoxels();
    irtkRealPixel *ptr2tgt  = target.GetPointerToVoxels();
    noOfPoints = target.GetNumberOfVoxels();

    for (m = 0; m < noOfPoints; m++){
      if (*ptr2tgt > Tp)
        *ptr2mask = 1;
      else
        *ptr2mask = 0;

      ++ptr2tgt;
      ++ptr2mask;
    }

  } else {
    mask.Read(mask_name);
  }

  // Make an identity global transformation.
  irtkAffineTransformation *trAffine = new irtkAffineTransformation;

  xdim = target.GetX();
  ydim = target.GetY();
  zdim = target.GetZ();

  if (regressAffine == true){
    // Estimate the global affine transformation.

    irtkPointSet targetPts;
    irtkPointSet sourcePts;

    // Collect point data.

    noOfPoints = target.GetNumberOfVoxels();

    int incr;
    incr = 1;
    while ((noOfPoints / incr) > MAX_PTS_PAREG){
      incr++;
    }

    // Subsample uniformly by increments of 'incr'
    noOfPoints = 0;
    ptID = -1;

    // Loop over all voxels.
    for (k = 0; k < zdim; k++){
      for (j = 0; j < ydim; j++){
        for (i = 0; i < xdim; i++){

          ptID++;
          // Should we sample it or not?
          if ((ptID % incr) != 0)
            continue;

          // Get two copies of current image coordinates.
          irtkPoint p(i, j, k);
          irtkPoint q(i, j, k);

          // Transform points into target world coordinates.
          target.ImageToWorld(p);
          target.ImageToWorld(q);

          targetPts.Add(p);

          // Transform one point to source coordinates.
          for (m = 0; m < noOfDofs; m++){
            transformation[m]->Transform(q);
          }

          sourcePts.Add(q);

          noOfPoints++;
        }
      }
    }


    // Estimate global affine component

    irtkPointAffineRegistration *pareg = new irtkPointAffineRegistration;
    // Set input and output for the registration filter
    irtkPointSet tmp1 = targetPts;
    irtkPointSet tmp2 = sourcePts;
    pareg->SetInput(&tmp1, &tmp2);
    pareg->SetOutput(trAffine);

    // Run registration filter
    pareg->Run();


  }

  // Finally can calculate distances.
  noOfPoints = 0;
  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        if (mask(i,j,k) > 0){
          noOfPoints++;
        }
      }
    }
  }

  irtkRealImage xdisp(target_name);
  irtkRealImage ydisp(target_name);
  irtkRealImage zdisp(target_name);

  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        if (mask(i,j,k) > 0){

          irtkPoint p(i, j, k);
          irtkPoint q(i, j, k);

          // Transform points into target world coordinates
          target.ImageToWorld(p);
          target.ImageToWorld(q);

          // Apply global Affine to one copy of the target
          // points (this is the identity if no regression was done).
          trAffine->irtkTransformation::Transform(p);

          // Transform the other copy by the input dofs.
          for (m = 0; m < noOfDofs; m++){
            transformation[m]->Transform(q);
          }

          xdisp(i, j, k) = p._x - q._x;
          ydisp(i, j, k) = p._y - q._y;
          zdisp(i, j, k) = p._z - q._z;
        }

      }
    }
  }


  vals  = new float[1 + noOfPoints];
  vals2 = new float[1 + noOfPoints];
  L  = new float[1 + noOfPoints];
  L2 = new float[1 + noOfPoints];

  sumD = sumD2 = sumL = sumL2 = 0.0;

  noOfPoints = 0;
  for (k = 1; k < zdim-1; k++){
    for (j = 1; j < ydim-1; j++){
      for (i = 1; i < xdim-1; i++){

        if (mask(i,j,k) > 0){

          // Calculate distance.
          x = xdisp(i, j, k);
          y = ydisp(i, j, k);
          z = zdisp(i, j, k);

          val = x*x + y*y + z*z;

          sumD2 += val;
          vals2[1 + noOfPoints] = val;

          val = sqrt(val);

          sumD += val;
          vals[1 + noOfPoints]  = val;

          dx = xdisp(i+1, j, k) - xdisp(i-1, j, k);
          dy = ydisp(i, j+1, k) - ydisp(i, j-1, k);
          dz = zdisp(i, j, k+1) - zdisp(i, j, k-1);

          val = x*x + y*y + z*z + alpha*alpha*(dx*dx + dy*dy + dz*dz);

          sumL2 += val;
          L2[1 + noOfPoints] = val;

          val = sqrt(val);

          sumL += val;
          L[1 + noOfPoints] = val;


          noOfPoints++;
        }

      }
    }
  }



  sort(noOfPoints, vals);
  sort(noOfPoints, vals2);
  sort(noOfPoints, L);
  sort(noOfPoints, L2);

  m = 1 + (int) round(0.5 * (noOfPoints - 1));

  // Means.
  cout << sumD / (double (noOfPoints));
  cout << " " << sumD2 / (double (noOfPoints)) << " ";
  // Medians and point count.
  cout << vals[m] << " " << vals2[m] << " " << noOfPoints << " ";

  cout << sumL / (double (noOfPoints));
  cout << " " << sumL2 / (double (noOfPoints)) << " ";
  cout << L[m] << " " << L2[m];
  cout << endl;

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
