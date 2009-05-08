#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "\timageabs [input] [output]" << endl;
  cerr << "\tFind absolute value of the input." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  irtkRealImage input;

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  input.Read(input_name);

  irtkRealPixel *ptr2pix = input.GetPointerToVoxels();
  int voxels = input.GetNumberOfVoxels();

  for (int i = 0; i < voxels; ++i){
    *ptr2pix = fabs(*ptr2pix);
    ++ptr2pix;
  }

  input.Write(output_name);

}



