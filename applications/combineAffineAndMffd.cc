#include <irtkTransformation.h>

char *affine_name = NULL;
char *mffd_name = NULL;
char *affine_out_name = NULL;

void usage()
{
  cerr << "Usage: combineAffineAndMffd [areg_in] [mffd_in] [mffd_out]" << endl;
  cerr << "Combine the [areg_in] with [mffd_in].  Control point values in mffd_in are" << endl;
  cerr << "multiplied by the jacobian of areg_in.  Write result to [mffd_out]."  << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 2)
    usage();

  int i, j, k, xdim, ydim, zdim;
  double x, y, z, xnew, ynew, znew;
  irtkMatrix jac(3, 3);
  irtkMatrix globalMatrix(4, 4);

  affine_name  = argv[1];
  argc--;
  argv++;
  mffd_name  = argv[1];
  argc--;
  argv++;
  affine_out_name  = argv[1];
  argc--;
  argv++;

  irtkAffineTransformation *affine;
  irtkMultiLevelFreeFormTransformation *mffd;
  affine = new irtkAffineTransformation;
  affine->irtkTransformation::Read(affine_name);
  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(mffd_name);

  if (mffd->NumberOfLevels() > 1){
    cerr << "Only implemented for a single level mffd." << endl;
    exit(1);
  }

  // The ffd that will need adjusting.
  irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));

  // Get the global affine matrix.
  globalMatrix = affine->GetMatrix();

  // Copy the upper left 3x3 block which is assumed to be the jacobian.
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      jac(i, j) = globalMatrix(i, j);
    }
  }
  // Assign the global part.
  mffd->irtkAffineTransformation::PutMatrix(globalMatrix);

  // Now adjust the control point values.
  xdim = ffd->GetX();
  ydim = ffd->GetY();
  zdim = ffd->GetZ();

  // Make the jacobian go from target to source.
  jac.Invert();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        // Get the control point values.
        ffd->Get(i, j, k, x, y, z);
        // Transform by the jacobian matrix.
        xnew = jac(0, 0) * x + jac(0, 1) * y + jac(0, 2) * z;
        ynew = jac(1, 0) * x + jac(1, 1) * y + jac(1, 2) * z;
        znew = jac(2, 0) * x + jac(2, 1) * y + jac(2, 2) * z;

        ffd->Put(i, j, k, xnew, ynew, znew);
      }
    }
  }

  mffd->irtkTransformation::Write(affine_out_name);

}


//   bool ok;

//   while (argc > 1){
//     ok = false;
//     if ((ok == false) && (strcmp(argv[1], "-out") == 0)){
//       argc--;
//       argv++;
//       dofout_name  = argv[1];
//       argc--;
//       argv++;
//       ok = true;
//     }
//     if (ok == false){
//       cerr << "Can not parse argument " << argv[1] << endl;
//       exit(1);
//     }
//   }

//   cout << "Before print" << endl;
//   cout << "Affine input" << endl;
//   affine->Print();
//   cout << "mffd input" << endl;
//   mffd->Print();

//   cout << "final mffd output print" << endl;
//   mffd->Print();
//   if (dofout_name != NULL){
//     mffd->Invert();
//     mffd->Write(dofout_name);
//     cout << "final output print" << endl;
//     mffd->Print();
//   }





