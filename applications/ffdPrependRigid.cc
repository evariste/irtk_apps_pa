
#include <irtkGeometry.h>
#include <irtkTransformation.h>

// Default filenames
char *rigid_in_name  = NULL;
char *ffd_in_name  = NULL;
char *ffd_out_name  = NULL;

void usage()
{
  cerr << " Usage: thisApp rigidIn mffdIn mffdOut" << endl;
  cerr << "  " << endl;
  cerr << " Input set up: " << endl;
  cerr << "  imageA ---------> ImageB --------> ImageC" << endl;
  cerr << "          rigidIn           mffdIn" << endl;
  cerr << "  " << endl;
  cerr << " I.e. imageB is the target for mffdIn and the source for rigidIn. " << endl;
  cerr << "  " << endl;
  cerr << " Output: " << endl;
  cerr << "  " << endl;
  cerr << "  imageA -------------------------> ImageC " << endl;
  cerr << "                  mffdOut" << endl;
  cerr << "  I.e. the rigid transformation is incorporated into the input MFFD to give the output MFFD." << endl;
  cerr << "  " << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int ok;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse source and target images
  rigid_in_name = argv[1];
  argc--;
  argv++;
  ffd_in_name = argv[1];
  argc--;
  argv++;
  ffd_out_name = argv[1];
  argc--;
  argv++;



  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-something") == 0)) {
      argc--;
      argv++;
//      some_variable = argv[1];
//      // or
//      some_variable = atoi(argv[1]);
//      // or
//      some_variable = atof(argv[1]);
//      // And possibly
//      argc--;
//      argv++;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  // imageA -------------> imageB ---------------> imageC
  //         x -> Rx = y           y-> Ay + u(y)
  // Or
  //         x -> Rx              Rx-> ARx + u(Rx)
  // I.e.
  //            x -> ARx + u(Rx)
  // In imageA terms:
  //            x -> Bx + uu(x)
  // So
  //       B = AR  , uu(x) = u(Rx)

  // R is a rigid transformation matrix.
  // A is an affine global part of the MFFD and u(x) is the local displacement field represented by a FFD.


  // Read input rigid.
  irtkRigidTransformation *rigidTransf = new irtkRigidTransformation;

  rigidTransf->irtkTransformation::Read(rigid_in_name);

  // Read input MFFD
  irtkMultiLevelFreeFormTransformation *mffd_in = new irtkMultiLevelFreeFormTransformation;
  mffd_in->irtkTransformation::Read(ffd_in_name);


  irtkBSplineFreeFormTransformation3D *ffd_in = dynamic_cast<irtkBSplineFreeFormTransformation3D *>(mffd_in->GetLocalTransformation(0));

  // Should probably check that mffd_in does not have more than one ffd in it!

  irtkMatrix mffd_in_global, rigidMat, invRigidMat, mffd_out_global;

  // Matrix A
  mffd_in_global = mffd_in->irtkAffineTransformation::GetMatrix();

  // Matrix R
  rigidMat = rigidTransf->GetMatrix();

  invRigidMat = rigidMat;
  invRigidMat.Invert();

  // A*R
  mffd_out_global = mffd_in_global * rigidMat;

  double x1, x2, y1, y2, z1, z2;
  double xaxis[3], yaxis[3], zaxis[3];
  double dx, dy, dz;

  ffd_in->GetSpacing(dx, dy, dz);

  ffd_in->GetOrientation(xaxis, yaxis, zaxis);

  // Get bounding box of input FFD
  x1 = 0.0;
  y1 = 0.0;
  z1 = 0.0;
  x2 = ffd_in->GetX() - 1;
  y2 = ffd_in->GetY() - 1;
  z2 = ffd_in->GetZ() - 1;

  ffd_in->LatticeToWorld(x1, y1, z1);
  ffd_in->LatticeToWorld(x2, y2, z2);

  // Transform to get the bounding box of the output FFD
  irtkVector v_in(4);
  irtkVector v_out(4);

  v_in(0) = x1;
  v_in(1) = y1;
  v_in(2) = z1;
  v_in(3) = 1;

  v_out = invRigidMat * v_in;

  x1 = v_out(0);
  y1 = v_out(1);
  z1 = v_out(2);

  v_in(0) = x2;
  v_in(1) = y2;
  v_in(2) = z2;
  v_in(3) = 1;

  v_out = invRigidMat * v_in;

  x2 = v_out(0);
  y2 = v_out(1);
  z2 = v_out(2);

  // Transform the the input axes to get the axes of the output FFD. Only need to rotate.

  irtkMatrix rotMat = invRigidMat;
  // Zero out the translation parameters.
  rotMat(0, 3) = 0;
  rotMat(1, 3) = 0;
  rotMat(2, 3) = 0;

  // Transform each axis in turn.

  v_in(0) = xaxis[0];
  v_in(1) = xaxis[1];
  v_in(2) = xaxis[2];

  v_out = rotMat * v_in;

  xaxis[0] = v_out(0);
  xaxis[1] = v_out(1);
  xaxis[2] = v_out(2);

  v_in(0) = yaxis[0];
  v_in(1) = yaxis[1];
  v_in(2) = yaxis[2];

  v_out = rotMat * v_in;

  yaxis[0] = v_out(0);
  yaxis[1] = v_out(1);
  yaxis[2] = v_out(2);

  v_in(0) = zaxis[0];
  v_in(1) = zaxis[1];
  v_in(2) = zaxis[2];

  v_out = rotMat * v_in;

  zaxis[0] = v_out(0);
  zaxis[1] = v_out(1);
  zaxis[2] = v_out(2);

  // Create the output FFD with the transformed control point lattice.
  irtkBSplineFreeFormTransformation3D *ffd_out = new irtkBSplineFreeFormTransformation3D(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);

  // Copy the control point parameters to the new FFD.
  int noOfDofs;
  int i;

  noOfDofs = ffd_in->NumberOfDOFs();
  for (i = 0; i < noOfDofs; ++i){
    ffd_out->Put(i, ffd_in->Get(i));
  }

  // Put the output FFD and global transformation into a new Multi Level FFD object.
  irtkMultiLevelFreeFormTransformation *mffd_out = new irtkMultiLevelFreeFormTransformation;

  mffd_out->PushLocalTransformation(ffd_out);
  mffd_out->PutMatrix(mffd_out_global);

  mffd_out->irtkTransformation::Write(ffd_out_name);


}
