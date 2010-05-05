#include <irtkImage.h>
#include <irtkTransformation.h>

// Input transformation
char *dofin_name = NULL;
char *affine_in_name = NULL;

// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffdputaffine [mffd_in] [areg_in] [mffd_out]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  // Check command line
  if (argc != 4){
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  affine_in_name = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transformation no. 1
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);

  irtkAffineTransformation *affineDof = new irtkAffineTransformation;
  affineDof->irtkTransformation::Read(affine_in_name);
  cout << "done" << endl;


  // `de-affine' the original ffd.
  int xdim, ydim, zdim, i, j, k;
  double x, y, z;
  double xAffCorr, yAffCorr, zAffCorr;

  xdim = mffd->GetLocalTransformation(0)->GetX();
  ydim = mffd->GetLocalTransformation(0)->GetY();
  zdim = mffd->GetLocalTransformation(0)->GetZ();

  // Correct the control point values in each of the input FFDs with
  // respect to their corresponding affine components.


  // Correcting the cp values requires the affine part to be S->T.
  irtkMatrix globalIn = mffd->irtkAffineTransformation::GetMatrix();
  globalIn.Invert();
  irtkFreeFormTransformation3D *ffd =
      	dynamic_cast<irtkFreeFormTransformation3D *> (mffd->GetLocalTransformation(0));
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        ffd->Get(i, j, k, x, y, z);

        xAffCorr = globalIn(0, 0) * x + globalIn(0, 1) * y + globalIn(0, 2) * z;
        yAffCorr = globalIn(1, 0) * x + globalIn(1, 1) * y + globalIn(1, 2) * z;
        zAffCorr = globalIn(2, 0) * x + globalIn(2, 1) * y + globalIn(2, 2) * z;

        ffd->Put(i, j, k, xAffCorr, yAffCorr, zAffCorr);
      }
    }
  }

  // apply the new affine matrix to the ffd.
  irtkMatrix globalOut = affineDof->GetMatrix();
  mffd->PutMatrix(globalOut);

  ffd =	dynamic_cast<irtkFreeFormTransformation3D *> (mffd->GetLocalTransformation(0));

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        ffd->Get(i, j, k, x, y, z);

        xAffCorr = globalOut(0, 0) * x + globalOut(0, 1) * y + globalOut(0, 2) * z;
        yAffCorr = globalOut(1, 0) * x + globalOut(1, 1) * y + globalOut(1, 2) * z;
        zAffCorr = globalOut(2, 0) * x + globalOut(2, 1) * y + globalOut(2, 2) * z;

        ffd->Put(i, j, k, xAffCorr, yAffCorr, zAffCorr);
      }
    }
  }


  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}


