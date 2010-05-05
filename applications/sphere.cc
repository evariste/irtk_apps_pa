////////////////////////////////////////////////////////////////
// DRAW A SPHERE.

#include <irtkImage.h>

char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << "sphere [input template image] [output image] [radius]" << endl;
  exit(1);
}

int main(int argc, char **argv){

  int radius;

  if (argc < 4){
    usage();
  }

  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;
  radius = atoi(argv[1]);
  argv++;
  argc--;

  irtkGreyImage img(input_name);

  double d;

  int xdim = img.GetX();
  int ydim = img.GetY();
  int zdim = img.GetZ();
  double centX = -0.5 + xdim / 2.0;
  double centY = -0.5 + ydim / 2.0;
  double centZ = -0.5 + zdim / 2.0;

  if (zdim > 1){
    for (int k = 0; k < zdim; ++k){
      for (int j = 0; j < ydim; ++j){
        for (int i = 0; i < xdim; ++i){
          d = round ( sqrt((double)((i - centX)*(i - centX) + (j - centY)*(j - centY) + (k - centZ)*(k - centZ))) ); 
          if (d <= radius){
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
        d = round ( sqrt((double)((i - centX)*(i - centX) + (j - centY)*(j - centY))) );
        if (d <= radius){
          img.Put(i, j, 0, 1);
        }else{
          img.Put(i, j, 0, 0);
        }
      }
    }
  }

  img.Write(output_name);
}
