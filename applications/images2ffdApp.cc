#include <irtkImage.h>

#include <irtkTransformation.h>

char *dofin_name, *dofout_name, *imagex_name, *imagey_name , *imagez_name;

void usage()
{
  cerr << "Usage: approximate [imagex] [imagey] [imagez] [dofin] [dofout] <options>\n";
  cerr << "where <options> is one or more of the following:\n" << endl;
  exit(1);
}

void resetFFD(irtkFreeFormTransformation *ffd){
  int i;
  int count = ffd->NumberOfDOFs();

  for (i = 0; i < count; ++i){
    ffd->PutStatus(i, _Active);
    ffd->Put(i, 0.0);
  }
}

void images2ffd(irtkRealImage *imagex, irtkRealImage *imagey, irtkRealImage *imagez,
                irtkFreeFormTransformation3D *ffd){

  int i, j, k;
  double px, py, pz;

  int voxels = imagex->GetNumberOfVoxels();

  int xdim = imagex->GetX();
  int ydim = imagex->GetY();
  int zdim = imagex->GetZ();

  // Allocate some memory
  double *x1 = new double[voxels];
  double *y1 = new double[voxels];
  double *z1 = new double[voxels];
  double *x2 = new double[voxels];
  double *y2 = new double[voxels];
  double *z2 = new double[voxels];

  int count = 0;

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){

        px = i;
        py = j;
        pz = k;

        imagex->ImageToWorld(px, py, pz);

        x1[count] = px;
        y1[count] = py;
        z1[count] = pz;

        x2[count] = imagex->Get(i, j, k);
        y2[count] = imagey->Get(i, j, k);
        z2[count] = imagez->Get(i, j, k);

        ++count;

      }
    }
  }

  // Approximate transformations
  ffd->Approximate(x1, y1, z1, x2, y2, z2, voxels);

  delete [] x1;
  delete [] y1;
  delete [] z1;
  delete [] x2;
  delete [] y2;
  delete [] z2;

}

int main(int argc, char **argv)
{
  bool ok;

  // Check command line
  if (argc < 6){
    usage();
  }

  imagex_name = argv[1];
  argv++;
  argc--;
  imagey_name = argv[1];
  argv++;
  argc--;
  imagez_name = argv[1];
  argv++;
  argc--;
  dofin_name  = argv[1];
  argv++;
  argc--;
  dofout_name = argv[1];
  argv++;
  argc--;

  // Read images.
  irtkRealImage *imagex = new irtkRealImage(imagex_name);
  irtkRealImage *imagey = new irtkRealImage(imagey_name);
  irtkRealImage *imagez = new irtkRealImage(imagez_name);

  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  // Only need the ffd for its control point locations.
  irtkFreeFormTransformation3D *ffd =
	  dynamic_cast<irtkFreeFormTransformation3D *>
		  (mffd->GetLocalTransformation(0));


  irtkMultiLevelFreeFormTransformation *mffdOut = new irtkMultiLevelFreeFormTransformation;

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  resetFFD(ffd);
  images2ffd(imagex, imagey, imagez, ffd);

  // Write transformation collection
  mffdOut->PushLocalTransformation(ffd);
  mffdOut->irtkTransformation::Write(dofout_name);

}
