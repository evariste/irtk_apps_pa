#include <irtkImage.h>

#define MAX 1000

int main(int argc, char **argv)
{
  int x, y, z, i, n;
  irtkGreyPixel label[MAX];
  irtkGreyImage image;

  if (argc < 4){
    cerr << "Usage: extractlabel [input] [output] [label]" << endl;
    exit(1);
  }

  // Read image
  image.Read(argv[1]);

  // Read label
  n = 0;
  for (i = 0; i < argc-3; i++){
    label[i] = atoi(argv[i+3]);
    n++;
  }

  // Extract label
  irtkGreyImage image2(image);
  for (z = 0; z < image.GetZ(); z++){
    for (y = 0; y < image.GetY(); y++){
      for (x = 0; x < image.GetX(); x++){
	image2(x, y, z) = 0;
      }
    }
  }

  for (z = 0; z < image.GetZ(); z++){
    for (y = 0; y < image.GetY(); y++){
      for (x = 0; x < image.GetX(); x++){
	for (i = 0; i < n; i++){
	  if (image(x, y, z) == label[i]){
	    image2(x, y, z) = 1;
	    break;
	  }
	}
      }
    }
  }

  // Write image
  image2.Write(argv[2]);

}
