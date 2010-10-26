#include <irtkImage.h>
#include <irtkGaussianBlurring.h>
#include <nr.h>

char *input_name  = NULL;
char *output_name = NULL;

void usage()
{
  cerr << "Usage: equalise [input] [output] <-pad> <-sigma> <-max>" << endl;
  cerr << "output is a version of input with equalised histogram." << endl;
  cerr << "Use pad to focus on a region." << endl;
  cerr << "Output is scaled to range 0 to max (default original max value in image)." << endl;
  cerr << "Specify sigma if blurring required." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  irtkGreyImage mask;
  irtkRealImage input;

  int pad = MIN_GREY;
  double sigma = 0;

  irtkGreyPixel *pMask;
  irtkRealPixel *pIn;

  int voxels, i, j, unpaddedCount = 0;
  bool ok;
  float *intensities, *offsets, *ranks;
  irtkRealPixel maxVal, minVal;
  int offset, intervalStart, intervalEnd;

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  input.Read(input_name);
  mask.Read(input_name);

  input.GetMinMax(&minVal, &maxVal);

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;      argv++;
      pad = atoi(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;      argv++;
      sigma = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-max") == 0)){
      argc--;      argv++;
      maxVal = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  if (sigma > 0){
    // Blur image.
    irtkGaussianBlurringWithPadding<irtkRealPixel> gaussianBlurring(sigma, pad);
    gaussianBlurring.SetInput (&input);
    gaussianBlurring.SetOutput(&input);
    gaussianBlurring.Run();
  }

  // Count valid voxels.
  pMask  = mask.GetPointerToVoxels();
  voxels = mask.GetNumberOfVoxels();
  unpaddedCount = 0;

  for (i = 0; i < voxels; ++i){
    if (*pMask > pad){
      ++unpaddedCount;
    }
    ++pMask;
  }

  // Put all valid voxel intensities into an array.
  // Numerical recipes arrays are 1-indexed, so need an extra element at
  // the end.
  intensities = new float[unpaddedCount + 1];
  offsets     = new float[unpaddedCount + 1];
  ranks       = new float[unpaddedCount + 1];

  pMask = mask.GetPointerToVoxels();
  pIn   = input.GetPointerToVoxels();
  unpaddedCount = 0;
  offset = 0;

  for (i = 0; i < voxels; ++i){
    if (*pMask > pad){
      intensities[unpaddedCount] = *pIn;
      offsets[unpaddedCount]     = offset;
      ++unpaddedCount;
    }
    ++pMask;
    ++pIn;
    offset++;
  }

  // Prepare for nr 1-indexing
  for (i = unpaddedCount; i > 0; i--){
    intensities[i] = intensities[i - 1];
    offsets[i] = offsets[i - 1];
  }
  // This call sorts the intensities array, and re-orders indices in the
  // same way.
  sort2(unpaddedCount, intensities, offsets);
  // Undo the 1-indexing shift
  for(i = 0; i < unpaddedCount; i++){
    intensities[i] = intensities[i + 1];
    offsets[i] = offsets[i + 1];
  }

  // Work out the rank of each voxel from the shuffled indices.
  intervalStart =  0;
  for (i = 1; i < unpaddedCount; ++i){

    if (round(intensities[i]) > round(intensities[i - 1])){
      intervalEnd = i - 1;
      for (j = intervalStart; j <= intervalEnd; j++){
        ranks[j] = 0.5 * (intervalStart + intervalEnd);
      }
      intervalStart = i;
    }
  }
  if (intervalEnd < unpaddedCount - 1){
    for (j = intervalStart; j < unpaddedCount; j++)
      ranks[j] = 0.5 * (intervalStart + unpaddedCount - 1);
  }

  // Assign rank-based values to the output image.
  pIn = input.GetPointerToVoxels();

  for (i = 0; i < unpaddedCount; ++i){
    offset = (int) offsets[i];
    *(pIn + offset) = round(maxVal * ranks[i] / ((double) unpaddedCount));
  }

//   for (i = 0; i < voxels; ++i){
//     if (*pIn <= pad){
//       *pIn = 0;
//     }
//     ++pIn;
//   }

  input.Write(output_name);
}

