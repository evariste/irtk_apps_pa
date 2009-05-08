///////////////////////////////////////////////////////////////
// Exp of an image.
#include <irtkImage.h>

char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << "imageexp [input] [output]" << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 3)
    usage();

  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  irtkRealImage *imgin  = new irtkRealImage(input_name);
  irtkRealImage *imgout = new irtkRealImage(input_name);

  irtkRealPixel *ptr2in, *ptr2out;
  double min = FLT_MAX;
  int i, voxels;

  voxels = imgin->GetNumberOfVoxels();
  ptr2in = imgin->GetPointerToVoxels();
  ptr2out = imgout->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){

    *ptr2out = exp(*ptr2in);

    ++ptr2in;
    ++ptr2out;
  }

//   ptr2in = imgin->GetPointerToVoxels();
//   ptr2out = imgout->GetPointerToVoxels();
//   min = floor(min);

//   for (i = 0; i < voxels; ++i){
//     if (*ptr2in <= 0){
//       *ptr2out = min;
//     }
//     ++ptr2in;
//     ++ptr2out;
//   }

  imgout->Write(output_name);

}
