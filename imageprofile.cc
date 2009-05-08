#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: reflect [in]  [-u val] [-v val]" << endl;
  cerr << "Where -u and -v are distinct and each are one of -x, -y, -z" << endl;
  cerr << "the values give the image coordinates of the profile." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkRealImage image;
  int i, j, k, ok;
  int px, py, pz;
  px = py = pz = -1;

  if (argc != 6){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      px = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      py = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-z") == 0)){
      argc--;
      argv++;
      pz = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }

    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  // Read input
  image.Read(input_name);

  // Checking.
  if (px == -1){
    if (py < 0 || py > image.GetY() - 1 || pz < 0 || pz > image.GetZ() - 1){
      cerr << "y / z values out of range." << endl;
      exit(1);
    }
  }

  if (py == -1){
    if (px < 0 || px > image.GetX() - 1 || pz < 0 || pz > image.GetZ() - 1){
      cerr << "x / z values out of range." << endl;
      exit(1);
    }
  }

  if (pz == -1){
    if (py < 0 || py > image.GetY() - 1 || px < 0 || px > image.GetX() - 1){
      cerr << "y / x values out of range." << endl;
      exit(1);
    }
  }

  // Doing.
  if (px == -1){
    for (i = 0; i < image.GetX(); ++i){
      cout << image.Get(i, py, pz) << endl;
    }
  }

  if (py == -1){
    for (j = 0; j < image.GetY(); ++j){
      cout << image.Get(px, j, pz) << endl;
    }
  }

  if (pz == -1){
    for (k = 0; k < image.GetZ(); ++k){
      cout << image.Get(px, py, k) << endl;
    }
  }

  return 0;
}
