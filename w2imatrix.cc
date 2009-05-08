#include <irtkImage.h>

char *input_name = NULL;

void usage()
{
  cerr << "Usage: w2imatrix [image] [matrix name]\n";
  exit(1);
}


int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  irtkRealImage image;
  image.Read(argv[1]);
  image.GetWorldToImageMatrix().Write(argv[2]);

}
