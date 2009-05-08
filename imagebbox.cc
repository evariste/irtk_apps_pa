#include <irtkImage.h>

char *input_name = NULL;

void usage()
{
  cerr << "\timagebbox [input]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  int xdim, ydim, zdim;
  double x, y, z;

  irtkGreyImage input;

  input_name  = argv[1];
  argc--;
  argv++;

  input.Read(input_name);

  xdim = input.GetX();
  ydim = input.GetY();
  zdim = input.GetZ();

  x = 0;
  y = 0;
  z = 0;

  input.ImageToWorld(x, y, z);
  cout << x << " " << y << " " << z << endl;

  x = xdim - 1;
  y = ydim - 1;
  z = zdim - 1;

  input.ImageToWorld(x, y, z);
  cout << x << " " << y << " " << z << endl;

}



