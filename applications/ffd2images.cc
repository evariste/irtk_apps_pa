#include <irtkImage.h>
#include <irtkTransformation.h>

// Input transformation
char *dofin_name     = NULL;
char *templateName   = NULL;
char *outputBaseName = NULL;

// Output transformation
char *image_out_nameX = NULL;
char *image_out_nameY = NULL;
char *image_out_nameZ = NULL;

void usage()
{
  cerr << "Usage: ffd2images [ffdIn] [templateImage] [imageOutBaseName] <-options>\n" << endl;
  cerr << "Writes out x, y and z displacements at the voxel locations of the template image" << endl;
  cerr << "as given by the ffd." << endl;
  cerr << "Options: " << endl;
  cerr << "-mag       : Magnitude only." << endl;
  cerr << "-t [val]   : Time value in the case of a 4D transformation (default 0)." << endl;
  cerr << "-pad [val] : Padding value in template image." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k;
  double x, y, z, t;
  bool ok;
  bool magnitudeOnly = false;
  int xdim, ydim, zdim;
  char buffer[256];
  double val;
  double pad = -1.0 * FLT_MAX;

  // Check command line
  if (argc < 4){
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  templateName = argv[1];
  argc--;
  argv++;
  outputBaseName = argv[1];
  argc--;
  argv++;

  t = 0;

  // Parse options.
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-mag") == 0)){
      argc--;
      argv++;
      magnitudeOnly = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;
      argv++;
      pad = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-t") == 0)){
      argc--;
      argv++;
      t = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  // Read transformation no. 1
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  cout << "done" << endl;

  irtkRealImage *templateImg = new irtkRealImage(templateName);

  irtkRealImage *dispX  = new irtkRealImage(templateName);
  irtkRealImage *dispY  = new irtkRealImage(templateName);
  irtkRealImage *dispZ  = new irtkRealImage(templateName);
  irtkRealImage *magImg = new irtkRealImage(templateName);

  xdim = templateImg->GetX();
  ydim = templateImg->GetY();
  zdim = templateImg->GetZ();

  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        if (templateImg->Get(i, j, k) <= pad){
          continue;
        }

        x = i;
        y = j;
        z = k;
        templateImg->ImageToWorld(x, y, z);

        mffd->LocalDisplacement(x, y, z, t);

        dispX->Put(i, j, k, x);
        dispY->Put(i, j, k, y);
        dispZ->Put(i, j, k, z);
        val = sqrt((x*x) + (y*y) + (z*z));
        magImg->Put(i, j, k, val);
      }
    }
  }

  if (magnitudeOnly){
    sprintf(buffer, "%s.nii.gz", outputBaseName);
    magImg->Write(buffer);
  } else {
    sprintf(buffer, "%s_X.nii.gz", outputBaseName);
    dispX->Write(buffer);
    sprintf(buffer, "%s_Y.nii.gz", outputBaseName);
    dispY->Write(buffer);
    sprintf(buffer, "%s_Z.nii.gz", outputBaseName);
    dispZ->Write(buffer);
  }

}
