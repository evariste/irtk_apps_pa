
#include <irtkImage.h>

#include <irtkTransformation.h>

// Default filenames
char *input_name = NULL, *output_name, *dof_name  = NULL;


void usage()
{
  cerr << "Usage: ffdbending [imageIn] [imageOut] [ffd]\n" << endl;
  cerr << "<-padding value>         Padding value" << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, m, ok;
  
  double x, y, z, val, padding;
  int noOfLevels;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Initialize padding value
  padding = -1.0 * FLT_MAX;

 
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkRealImage *image = new irtkRealImage(input_name);
  cout << "done" << endl;

  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd;
  mffd = new irtkMultiLevelFreeFormTransformation;
  
  mffd->irtkTransformation::Read(dof_name);
  
  noOfLevels = mffd->NumberOfLevels();

  m = 0;
  n = 0;
  val = 1;
 
  for (k = 0; k < image->GetZ(); k++) {
    for (j = 0; j < image->GetY(); j++) {
      for (i = 0; i < image->GetX(); i++) {

        if (image->Get(i, j, k) > padding) {
          x = i;
          y = j;
          z = k;
        
          image->ImageToWorld(x, y, z);
          val = mffd->Bending(x, y, z);
          
        }
        image->Put(i, j, k, val);
      }
    }
  }
 
  // Write the final transformation estimate
  image->Write(output_name);
}
