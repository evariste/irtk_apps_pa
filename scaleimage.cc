#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "\tscaleimage [input] [scale factor] [output]" << endl;
  cerr << "\tScale intensity of each voxel in [input] by" << endl;
  cerr << "\t[scale factor] (fractional allowed)." << endl;
  cerr << "\tResult written to output." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 4){
    usage();
  }

  irtkRealImage input;
  irtkRealImage output;
  double scaleFactor;

  input_name  = argv[1];
  argc--;
  argv++;
  scaleFactor = atof(argv[1]);
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  input.Read(input_name);
  output.Read(input_name);

  irtkRealPixel *pPix = input.GetPointerToVoxels();
  irtkRealPixel *pPixOut = output.GetPointerToVoxels();

  int voxels = input.GetNumberOfVoxels();

  for (int i = 0; i < voxels; ++i){
    *pPixOut = (*pPix) * scaleFactor;
    ++pPix;
    ++pPixOut;
  }

  output.Write(output_name);
}



