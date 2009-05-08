#include <irtkImage.h>

// Default filenames
char *output_name = NULL, *input_name = NULL;

void usage()
{
  cerr << "Usage: fill [input] [output] n" << endl;
  cerr << "Input assumed to be a binary mask where background is zero and less.\n";
}

int main(int argc, char **argv)
{

  int i, j, k, nVox, count;
  int xdim, ydim, zdim;
  int xoffset, yoffset, zoffset;

  irtkGreyImage input;
  irtkGreyImage output;
  irtkGreyPixel *pPix1, *pPix2;

  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse images
  input_name = argv[1];
  argc--;
  argv++;

  output_name = argv[1];
  argc--;
  argv++;

  input.Read(input_name);
  output.Read(input_name);

  // Set all output pixels to one.
  pPix1 = output.GetPointerToVoxels();
  nVox = output.GetNumberOfVoxels();
  for (i = 0; i < nVox; ++i){
    *pPix1 = 1;
    ++pPix1;
  }



  xdim = input.GetX();
  ydim = input.GetY();
  zdim = input.GetZ();

  xoffset = 1;
  yoffset = input.GetX();
  zoffset = input.GetX() * input.GetY();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      pPix1 = input.GetPointerToVoxels(0, j, k);
      pPix2 = output.GetPointerToVoxels(0, j, k);

      count = 0;
      while (*pPix1 <= 0 && count < xdim){
        *pPix2 = 0;
        pPix1 += xoffset;
        pPix2 += xoffset;
        ++count;
      }

      pPix1 = input.GetPointerToVoxels(xdim - 1, j, k);
      pPix2 = output.GetPointerToVoxels(xdim - 1, j, k);
      count = 0;
      while (*pPix1 <= 0 && count < xdim){
        *pPix2 = 0;
        pPix1 -= xoffset;
        pPix2 -= xoffset;
        ++count;
      }
    }
  }

  for (k = 0; k < zdim; ++k){
    for (i = 0; i < xdim; ++i){
      pPix1 = input.GetPointerToVoxels(i, 0, k);
      pPix2 = output.GetPointerToVoxels(i, 0, k);
      count = 0;
      while (*pPix1 <= 0 && count < ydim){
        *pPix2 = 0;
        pPix1 += yoffset;
        pPix2 += yoffset;
        ++count;
      }

      pPix1 = input.GetPointerToVoxels(i, ydim - 1, k);
      pPix2 = output.GetPointerToVoxels(i, ydim - 1, k);
      count = 0;
      while (*pPix1 <= 0 && count < ydim){
        *pPix2 = 0;
        pPix1 -= yoffset;
        pPix2 -= yoffset;
        ++count;
      }
    }
  }

  for (j = 0; j < ydim; ++j){
    for (i = 0; i < xdim; ++i){
      pPix1 = input.GetPointerToVoxels(i, j, 0);
      pPix2 = output.GetPointerToVoxels(i, j, 0);
      count = 0;
      while (*pPix1 <= 0 && count < zdim){
        *pPix2 = 0;
        pPix1 += zoffset;
        pPix2 += zoffset;
        ++count;
      }

      pPix1 = input.GetPointerToVoxels(i, j, zdim - 1);
      pPix2 = output.GetPointerToVoxels(i, j, zdim - 1);
      count = 0;
      while (*pPix1 <= 0 && count < zdim){
        *pPix2 = 0;
        pPix1 -= zoffset;
        pPix2 -= zoffset;
        ++count;
      }
    }
  }

  output.Write(output_name);

}

