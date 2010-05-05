#include <irtkImage.h>

char *input_name = NULL;

void usage()
{
  cerr << "Usage: i2wmatrix [image] [matrix name]\n";
  exit(1);
}


int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  irtkRealImage image;
  image.Read(argv[1]);
  image.GetImageToWorldMatrix().Write(argv[2]);

}
