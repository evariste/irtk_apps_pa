////////////////////////////////////////////////////////////////
// DRAW A CUBOID.

#include <irtkImage.h>

char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << "cuboid [input template image] [output image]  [x1 y1 z1 x2 y2 z2]" << endl;
  cerr << "fill a cuboid with corners  x1 y1 z1 x2 y2 z2" << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv){

  int x1, y1, z1, x2, y2, z2;

  if (argc < 9){
    usage();
  }

  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  x1  = atoi(argv[1]);
  argv++;
  argc--;
  y1 = atoi(argv[1]);
  argv++;
  argc--;
  z1 = atoi(argv[1]);
  argv++;
  argc--;
  x2 = atoi(argv[1]);
  argv++;
  argc--;
  y2 = atoi(argv[1]);
  argv++;
  argc--;
  z2 = atoi(argv[1]);
  argv++;
  argc--;

  irtkGreyImage img(input_name);

  int xdim = img.GetX();
  int ydim = img.GetY();
  int zdim = img.GetZ();

  if (x1 > x2)
    swap(x1, x2);
  if (y1 > y2)
    swap(y1, y2);
  if (z1 > z2)
    swap(z1, z2);

  x1--;
  y1--;
  z1--;
  x2++;
  y2++;
  z2++;

  if (zdim > 1){
    for (int k = 0; k < zdim; ++k){
      for (int j = 0; j < ydim; ++j){
        for (int i = 0; i < xdim; ++i){

          if (i > x1 && i < x2 && j > y1 && j < y2 && k > z1 && k < z2){
            img.Put(i, j, k, 1);
          }else{
            img.Put(i, j, k, 0);
          }
        }
      }
    }
  } else {
    for (int j = 0; j < ydim; ++j){
      for (int i = 0; i < xdim; ++i){
        if (i > x1 && i < x2 && j > y1 && j < y2){
          img.Put(i, j, 0, 1);
        }else{
          img.Put(i, j, 0, 0);
        }
      }
    }
  }

  img.Write(output_name);
}
