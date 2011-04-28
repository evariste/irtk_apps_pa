#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>

// Default filenames
char *dofOut_name = NULL;

char **dof_name  = NULL;

#define MAX_DOFS 50
#define MAX_PTS_PAREG 10000

void usage()
{
  cerr << " Usage: ffdcomposeN [dofOut] <options>\n" << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " " << endl;
  cerr << " <-dofin file>      Transformation file. Multiple transformations can be given by" << endl;
  cerr << "                    repeatedly using this flag. They are processed in order they passed." << endl;
  cerr << " <-dofin_i file>    A transformation whose inverse should be applied." << endl;
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

  int i, j, k, m, xdim, ydim, zdim, noOfDofs, numberOfCPs, count;
  double x, y, z;
  double xStore, yStore, zStore;

  bool ok;

  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse dofOut_name
  dofOut_name = argv[1];
  argc--;
  argv++;



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


    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

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
  mffd_out->irtkTransformation::Read(dof_name[0]);  // takes the first DOF

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
  for (i = 0; i < xdim; i++){
    for (j = 0; j < ydim; j++){
      for (k = 0; k < zdim; k++){

        affd_out->Put(i, j, k, 0, 0, 0);

        x = i; // why? i not used later
        y = j; // why? j not used later
        z = k; // why? k not used later

        // Transform point from lattice coordinates to target coordinates
        affd_out->LatticeToWorld(x, y, z);

        xStore = x;
        yStore = y;
        zStore = z;

        // Transform point
        for (m = 0; m < noOfDofs; m++){
          if (invert[m])
            transformation[m]->Inverse(x, y, z);
          else
            transformation[m]->Transform(x, y, z);
        }

        // New location of the control point.
        xdata[count] = x - xStore;
        ydata[count] = y - yStore;
        zdata[count] = z - zStore;

        ++count; 

      }
    }
  }


  // Estimate an affine component for the collected displacements.

  // Make an identity global transformation.
  irtkAffineTransformation *trAffine = new irtkAffineTransformation;

  irtkPointSet targetPts;
  irtkPointSet sourcePts;

  // Collect point data.

  int noOfPoints;
  noOfPoints = affd_out->NumberOfDOFs() / 3;

  int incr = 1;
  while ((noOfPoints / incr) > MAX_PTS_PAREG){
    incr++;
  }

  // Subsample uniformly by increments of 'incr'

  count = -1;

  // Loop over all control points.
  for (i = 0; i < xdim; i++){
    for (j = 0; j < ydim; j++){
      for (k = 0; k < zdim; k++){

        count++;


        // Should we sample it or not?
        if ((count % incr) != 0)
          continue;

        // Get two copies of current image coordinates.
        x = i; y = j; z = k;
        irtkPoint p(x, y, z);
        irtkPoint q(x, y, z);

        // Transform points into target world coordinates.
        affd_out->LatticeToWorld(p);
        affd_out->LatticeToWorld(q);


        // The starting point
        targetPts.Add(p);

        // Where does p go to under all the compositions?
        q._x += xdata[count];
        q._y += ydata[count];
        q._z += zdata[count];

        sourcePts.Add(q);
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


  cout << "Estimated global affine component with following parameters:" << endl;
  trAffine->Print();


  // Remove the global affine part from the estimated displacements
  count = 0;
  for (i = 0; i < xdim; i++){
    for (j = 0; j < ydim; j++){
      for (k = 0; k < zdim; k++){

        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target coordinates
        affd_out->LatticeToWorld(x, y, z);

        xStore = x;
        yStore = y;
        zStore = z;

        trAffine->Transform(x, y, z);


        // New location of the control point.
        xdata[count] -= (x - xStore);
        ydata[count] -= (y - yStore);
        zdata[count] -= (z - zStore);

        ++count;

      }
    }
  }




  // Interpolate the ffd and write dof
  affd_out->Interpolate(xdata, ydata, zdata);


  mffd_out->PutMatrix(trAffine->GetMatrix());


  (void) mffd_out->PopLocalTransformation();
  mffd_out->PushLocalTransformation(affd_out);
  mffd_out->irtkTransformation::Write(dofOut_name);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;
  delete [] dof_name;
  delete [] invert;

}

