///////////////////////////////////////////////////////////////
// Centre of gravity of a mask image.

#include <irtkImage.h>

char *image_name = NULL;

void usage()
{
  cerr << "Usage: imagecentreofgravity  [image] <options>" << endl;
  cerr << " -pad : padding value , default = 0" << endl;
  exit(1);
}


int main(int argc, char **argv){

  if (argc < 2)
    usage();

  image_name = argv[1];
  argv++;
  argc--;

  irtkRealPixel *ptr2voxel;

  irtkRealImage *image = new irtkRealImage(image_name);

  int pad = 0;
  int ok = True;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-pad") == 0)){
      argc--;
      argv++;
      pad = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  double x, y, z;
  double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
  int i, j, k, count = 0;

  ptr2voxel = image->GetPointerToVoxels();

  for (k = 0; k < image->GetZ(); ++k){
    for (j = 0; j < image->GetY(); ++j){
      for (i = 0; i < image->GetX(); ++i){

        if (image->Get(i, j, k) > pad){
          x = i;
          y = j;
          z = k;
          image->ImageToWorld(x, y, z);

          sumX += x;
          sumY += y;
          sumZ += z;
          ++count;
        }


      }
    }

  }

  cout << sumX / (double (count)) << " ";
  cout << sumY / (double (count)) << " ";
  cout << sumZ / (double (count)) << endl;

}

