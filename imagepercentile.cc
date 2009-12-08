#include <irtkImage.h>

#include <nr.h>

#define MAXVALS 100

char *input_name = NULL;
char *mask_name = NULL;

void usage()
{
  cerr << " imagepercentile  [input] [N] [prctile_1] ... [prctile_N] <-options>" << endl;
  cerr << " Find percentiles of an image." << endl;
  cerr << " N : Number of percentiles required (max " << MAXVALS << ")." << endl;
  cerr << " [prctile_1] ... [prctile_N] : percentiles required." << endl;
  cerr << " Options:" << endl;
  cerr << " -pad : Ignored padding value." << endl;
  cerr << " -q   : Quiet; space separated percentiles only." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  irtkRealImage input;
  int i, numberOfPercentiles;
  int percentile[MAXVALS];
  int quiet = False;

  double pad = -1.0 * FLT_MAX;
  int ok;
  int voxels;
  int count, index;

  input_name  = argv[1];
  argc--;
  argv++;
  numberOfPercentiles = atoi(argv[1]);
  argc--;
  argv++;

  if (numberOfPercentiles > MAXVALS){
    cerr << "Cannot have more than " << MAXVALS << " percentiles!" << endl << endl;
    usage();
    exit(1);
  }

  for (i = 0; i < numberOfPercentiles; i++){
    percentile[i] = atoi(argv[1]);
    argc--;
    argv++;
  }

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-pad") == 0)){
      argc--;
      argv++;
      pad = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  input.Read(input_name);
  irtkRealPixel *pPix;
  irtkGreyPixel *pMask;

  voxels = input.GetNumberOfVoxels();

  irtkGreyImage mask;
  if (mask_name == NULL){
    mask.Read(input_name);

    pPix = input.GetPointerToVoxels();
    pMask = mask.GetPointerToVoxels();

    for (i = 0; i < voxels; ++i){
      if (*pPix > pad){
        *pMask = 1;
      } else {
        *pMask = 0;
      }
      ++pPix;
      ++pMask;
    }
  } else {
    mask.Read(mask_name);
    if ((mask.GetX() != input.GetX()) ||
        (mask.GetY() != input.GetY()) ||
        (mask.GetZ() != input.GetZ())){
      cerr << "imagepercentile: Input and mask dimensions mismatch. Exiting." << endl;
      exit(1);
    }
  }
  
  pMask = mask.GetPointerToVoxels();
  
  // Find number of unpadded values.
  count = 0;
  for (i = 0; i < voxels; ++i){
    if (*pMask > 0){
      ++count;
    }
    ++pMask;
  }

  // NR subroutines use 1-indexing.
  float *data = new float[1 + count];
  
  pMask = mask.GetPointerToVoxels();
  pPix = input.GetPointerToVoxels();
  
  count = 0;
  
  for (i = 0; i < voxels; ++i){
    if (*pMask > 0){
      data[1 + count] = *pPix;
      ++count;
    }
    ++pPix;
    ++pMask;
  }

  sort(count, data);

  if (quiet == True){
    for (i = 0; i < numberOfPercentiles; i++){
      index = 1 + (int) round( (double) percentile[i] * (count - 1) / 100.0);
      cout << data[index] << " ";
    }
  } else {
    for (i = 0; i < numberOfPercentiles; i++){
      index = 1 + (int) round( (double) percentile[i] * (count - 1) / 100.0);
      cout << "  Percentile " << percentile[i] << " = " << data[index] << endl;
    }
  }

}
