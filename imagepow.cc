///////////////////////////////////////////////////////////////
// Exp of an image.
#include <irtkImage.h>

char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << "imagepow [input] [power] [output]" << endl;
  exit(1);
}

int main(int argc, char **argv){

  double power;

  if (argc < 3)
    usage();

  input_name = argv[1];
  argv++;
  argc--;
  power = atof(argv[1]);
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  irtkRealImage *imgin  = new irtkRealImage(input_name);
  irtkRealImage *imgout = new irtkRealImage(input_name);

  irtkRealPixel *ptr2in, *ptr2out;
  int i, voxels;

  voxels = imgin->GetNumberOfVoxels();
  ptr2in = imgin->GetPointerToVoxels();
  ptr2out = imgout->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){

    *ptr2out = pow(((double) *ptr2in), power);

    ++ptr2in;
    ++ptr2out;
  }


  imgout->Write(output_name);

}
