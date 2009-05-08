#include <irtkImage.h>

char *input_1 = NULL, *input_2 = NULL, *output_name = NULL;

void usage()
{
  cerr << "\taddimages [input 1] [input 2] [output]" << endl;
  cerr << "\tAdd images on a voxel by voxel basis." << endl;
  cerr << "\tImages must of course have the same dimensions." << endl;
  cerr << "\tResult written to output." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 4){
    usage();
  }

  irtkRealImage input1;
  irtkRealImage input2;
  irtkRealImage output;

  input_1  = argv[1];
  argc--;
  argv++;
  input_2  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  input1.Read(input_1);
  input2.Read(input_2);
  output.Read(input_1);

  int x1 = input1.GetX();
  int x2 = input2.GetX();
  int y1 = input1.GetY();
  int y2 = input2.GetY();
  int z1 = input1.GetZ();
  int z2 = input2.GetZ();

  if (x1 != x2 || y1 != y2 || z1 != z2){
    usage();
  }

  irtkRealPixel *pPix1 = input1.GetPointerToVoxels();
  irtkRealPixel *pPix2 = input2.GetPointerToVoxels();
  irtkRealPixel *pPixOut = output.GetPointerToVoxels();

  int voxels = input1.GetNumberOfVoxels();

  for (int i = 0; i < voxels; ++i){
    *pPixOut = *pPix1 + *pPix2;
    ++pPix1;
    ++pPix2;
    ++pPixOut;
  }

  output.Write(output_name);

}



