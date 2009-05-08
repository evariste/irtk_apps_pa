#include <irtkImage.h>
#include <irtkTransformation.h>

char *inputDofName       = NULL;
char *templateDofName     = NULL;
char *outputDofName       = NULL;

void usage(){
  cout << "\n Usage: ffdinvert [mffdIn] [templateMffd] [mffdOut]" << endl;
  cout << "Invert the FFD [mffdIn].  The result is estimated for the " << endl;
  cout << "control points defined in [templateMffd] and is written to [mffdOut]." << endl;
  cout << "For images A and B assume that: " << endl;
  cout << "  mffdIn       : A -> B  (i.e. B is the target)" << endl;
  cout << "  templateMffd : has a control point lattice defined using A as the target." << endl;
  cout << "  mffdOut      : B -> A" << endl;
  cout << "" << endl;
  exit(0);
}

int main(int argc, char **argv)
{
  irtkMatrix jac, global;
  int i, j, k, xdim, ydim, zdim, numberOfCPs, count;
  double x, y, z;

  if (argc != 4){
    usage();
  }

  inputDofName    = argv[1];
  argc--;
  argv++;
  templateDofName = argv[1];
  argc--;
  argv++;
  outputDofName   = argv[1];
  argc--;
  argv++;

  // Read first transformation, this transformation is to be
  // inverted.  mffd_in : A -> B
  irtkMultiLevelFreeFormTransformation *mffd_in = new irtkMultiLevelFreeFormTransformation;
  mffd_in->irtkTransformation::Read(inputDofName);

  // Read transformation that has the lattice where we want to
  // define the output.
  // mffd_template : B -> A
  irtkMultiLevelFreeFormTransformation *mffd_template = new irtkMultiLevelFreeFormTransformation;
  mffd_template->irtkTransformation::Read(templateDofName);

  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *affd_template =
	  dynamic_cast<irtkFreeFormTransformation3D *>(mffd_template->GetLocalTransformation(0));
  xdim = affd_template->GetX();
  ydim = affd_template->GetY();
  zdim = affd_template->GetZ();

  numberOfCPs = xdim * ydim * zdim;

  // Space to store the control point displacements.
  double *xdata = new double[numberOfCPs];
  double *ydata = new double[numberOfCPs];
  double *zdata = new double[numberOfCPs];

  // Store the global transformation.
  global = mffd_in->irtkAffineTransformation::GetMatrix();
  global.Invert();
  // global : B -> A, this will be the global part of the final
  // output mffd.

  count = 0;

  irtkMatrix cp(4, 1);
  irtkMatrix p(4, 1);

  // Loop over points in source lattice.
  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target
        // coordinates
        affd_template->LatticeToWorld(x, y, z);

        // Keep a copy of the control point location.
        cp(0, 0) = x;
        cp(1, 0) = y;
        cp(2, 0) = z;
        cp(3, 0) = 1;

        // How is the control point affected by the global
        // transformation B->A?
        p = global * cp;

        // What is the estimate for the full inverse
        // transformation (global and local)?
        mffd_in->Inverse(x, y, z);

        // Subtract the global transformation to get an estimate
        // for the displacement at the control point.
        xdata[count] = x - p(0, 0);
        ydata[count] = y - p(1, 0);
        zdata[count] = z - p(2, 0);

        ++count;
      }
    }
  }

  // Interpolate the ffd.
  affd_template->Interpolate(xdata, ydata, zdata);

  // The output DOF is now in mffd_template and its affine portion needs to be added.
  mffd_template->irtkAffineTransformation::PutMatrix(global);

  mffd_template->irtkTransformation::Write(outputDofName);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;

}

///////////////////////////////////////////////////////////////////////////////////////

// Function that is now unused
void FFDRemoveAffine(irtkMultiLevelFreeFormTransformation *mffd){

  // Correct the control points of an FFD with respect to the global
  // transformation of the containing MFFD.
  int i, j, k;
  int xdim, ydim, zdim;
  double x, y, z;
  double xAffCorr, yAffCorr, zAffCorr;
  irtkMatrix global(4, 4);

  if (mffd->NumberOfLevels() > 1){
    cerr << "The mffd to be transformed must have only one level." << endl;
    exit(1);
  }

  irtkFreeFormTransformation3D *affd =
	  dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));

  xdim = affd->GetX();
  ydim = affd->GetY();
  zdim = affd->GetZ();

  global = mffd->irtkAffineTransformation::GetMatrix();
  // Correcting the CP values requires the affine part to be S->T.
  global.Invert();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        affd->Get(i, j, k, x, y, z);

        xAffCorr = global(0, 0) * x + global(0, 1) * y + global(0, 2) * z;
        yAffCorr = global(1, 0) * x + global(1, 1) * y + global(1, 2) * z;
        zAffCorr = global(2, 0) * x + global(2, 1) * y + global(2, 2) * z;

        affd->Put(i, j, k, xAffCorr, yAffCorr, zAffCorr);

      }
    }
  }

}


