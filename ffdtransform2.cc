#include <irtkImage.h>
#include <irtkTransformation.h>

char *transformingDofName = NULL;
char *sourceDofName       = NULL;
char *templateDofName     = NULL;
char *outputDofName       = NULL;

void usage(){
  cout << "Usage: ffdtransform2 [transformingDof] [sourceDof] [templateDof] [outputDof]" << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;
  exit(0);
}

void FFDRemoveAffine(irtkMultiLevelFreeFormTransformation *);

int main(int argc, char **argv)
{
  irtkMatrix jac, global;
  int i, j, k, xdim, ydim, zdim, numberOfCPs, count;
  double x, y, z, v1[3];
  float v2[3];
  double val;

  if (argc != 5){
    usage();
  }

  transformingDofName = argv[1];
  argc--;
  argv++;
  sourceDofName       = argv[1];
  argc--;
  argv++;
  templateDofName     = argv[1];
  argc--;
  argv++;
  outputDofName       = argv[1];
  argc--;
  argv++;


  // Read first transformation, this transformation takes locations in the
  // target to locations in the source image space and is used to pull back
  // a transformation defined in the source image space.
  irtkMultiLevelFreeFormTransformation *mffd_transform = new irtkMultiLevelFreeFormTransformation;
  mffd_transform->irtkTransformation::Read(transformingDofName);

  // Read second transformation, this transformation is defined in
  // coordinate system of the source image for the previous transformation.
  irtkMultiLevelFreeFormTransformation *mffd_source = new irtkMultiLevelFreeFormTransformation;
  mffd_source->irtkTransformation::Read(sourceDofName);

  // Read transformation that has the lattice where we want to define the
  // output.
  irtkMultiLevelFreeFormTransformation *mffd_template = new irtkMultiLevelFreeFormTransformation;
  mffd_template->irtkTransformation::Read(templateDofName);

  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *affd_template =
	  dynamic_cast<irtkFreeFormTransformation3D *>
		   (mffd_template->GetLocalTransformation(0));

  xdim = affd_template->GetX();
  ydim = affd_template->GetY();
  zdim = affd_template->GetZ();

  numberOfCPs = xdim * ydim * zdim;

  // Space to store the control point displacements.
  double *xdata = new double[numberOfCPs];
  double *ydata = new double[numberOfCPs];
  double *zdata = new double[numberOfCPs];

  // Remove the effect of the global affine transformation from the
  // (single) FFD in mffd_source.
  FFDRemoveAffine(mffd_source);

  count = 0;

  // Loop over points in target
  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target coordinates
        affd_template->LatticeToWorld(x, y, z);

        // Calculate total jacobian at the point
        val = mffd_transform->irtkTransformation::Jacobian(x, y, z);

        if ((val > 0.1) && (val < 5)){

          mffd_transform->Jacobian(jac, x, y, z);
          jac.Invert();

          // Transform target world coordinates to source coordinates
          mffd_transform->Transform(x, y, z);

          // Store world coords
          v1[0] = x;
          v1[1] = y;
          v1[2] = z;

          // Tranform point (recall mffd_source has been affine corrected).
          mffd_source->LocalDisplacement(v1[0], v1[1], v1[2]);

          // Transform deformation vector
          v2[0] = jac(0,0) * v1[0] + jac(0,1) * v1[1] + jac(0,2) * v1[2];
          v2[1] = jac(1,0) * v1[0] + jac(1,1) * v1[1] + jac(1,2) * v1[2];
          v2[2] = jac(2,0) * v1[0] + jac(2,1) * v1[1] + jac(2,2) * v1[2];

        } else {
          v2[0] = 0;
          v2[1] = 0;
          v2[2] = 0;
        }

        xdata[count] = v2[0];
        ydata[count] = v2[1];
        zdata[count] = v2[2];
        ++count;
      }
    }
  }

  // Interpolate the ffd.
  affd_template->Interpolate(xdata, ydata, zdata);

  // The output DOF is now in mffd_template and its affine portion needs to be reset.
  global = mffd_template->irtkAffineTransformation::GetMatrix();
  global.Ident();
  mffd_template->irtkAffineTransformation::PutMatrix(global);

  mffd_template->irtkTransformation::Write(outputDofName);

}


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
	  dynamic_cast<irtkFreeFormTransformation3D *>
		  (mffd->GetLocalTransformation(0));

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

