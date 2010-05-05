#include <irtkImage.h>
#include <irtkTransformation.h>

// Input transformation
char *dofout_name = NULL;

// Output transformation
char *image_in_nameX = NULL;
char *image_in_nameY = NULL;
char *image_in_nameZ = NULL;

void usage()
{
  cerr << "Usage: ffd2images [imageInX] [imageInY] [imageInZ] [dofout] \n" << endl;
  exit(1);
}

bool sameSizeImages(irtkRealImage *a, irtkRealImage *b){
  if (a->GetX() != b->GetX() || a->GetY() != b->GetY() || a->GetZ() != b->GetZ() ){
    return false;
  }
  return true;
}

int main(int argc, char **argv)
{
  int i, j, k;
  double x, y, z;

  double x1, y1, z1, x2, y2, z2;
  double xsize, ysize, zsize;
  double xaxis[3], yaxis[3], zaxis[3];

  irtkRealImage *cpsX, *cpsY, *cpsZ;

  irtkMultiLevelFreeFormTransformation *mffd;
  irtkBSplineFreeFormTransformation *affd;

  // Check command line
  if (argc != 5){
    usage();
  }

  // Parse file names
  image_in_nameX  = argv[1];
  argc--;
  argv++;
  image_in_nameY  = argv[1];
  argc--;
  argv++;
  image_in_nameZ  = argv[1];
  argc--;
  argv++;

  dofout_name  = argv[1];
  argc--;
  argv++;

  // Read images.
  cpsX = new irtkRealImage(image_in_nameX);
  cpsY = new irtkRealImage(image_in_nameY);
  cpsZ = new irtkRealImage(image_in_nameZ);

  if (( ! sameSizeImages(cpsX , cpsY) ) || ( ! sameSizeImages(cpsX , cpsZ))){
    cerr << "Images must have same dimensions." << endl;
    exit(1);
  }

  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = cpsX->GetX() - 1;
  y2 = cpsX->GetY() - 1;
  z2 = cpsX->GetZ() - 1;

  cpsX->ImageToWorld(x1, y1, z1);
  cpsX->ImageToWorld(x2, y2, z2);

  cpsX->GetPixelSize(&xsize, &ysize, &zsize);
  cpsX->GetOrientation(xaxis, yaxis, zaxis);

  affd = new irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, xsize, ysize, zsize, xaxis, yaxis, zaxis);

  for (i = 0; i < cpsX->GetX(); i++){
    for (j = 0; j < cpsX->GetY(); j++){
      for (k = 0; k < cpsX->GetZ(); k++){

        x = cpsX->Get(i, j, k);
        y = cpsY->Get(i, j, k);
        z = cpsZ->Get(i, j, k);

        affd->Put(i, j, k, x, y, z);

      }
    }
  }

  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->PushLocalTransformation(affd);
  mffd->irtkTransformation::Write(dofout_name);

}
