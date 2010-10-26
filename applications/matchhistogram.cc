#include <irtkImage.h>
#include <irtkGaussianBlurring.h>
#include <nr.h>

char *ref_name    = NULL;
char *input_name  = NULL;
char *output_name = NULL;

void usage()
{
  cerr << "Usage: matchhistogram [reference] [input] [output] <options>" << endl;
  cerr << "[output] is a version of [input] with the " << endl;                    
  cerr << "histogram matched to the histogram of [reference]." << endl;                    
  cerr << "I.e. the cumulative distributions should be" << endl;                    
  cerr << " similar." << endl;
  cerr << "Options: " << endl;
  cerr << "-pad <value>   : Ignore values less than equal to" << endl;                    
  cerr << "                 the given value." << endl;
  cerr << "-sigma <value> : Blur images with Gaussian kernel " << endl;
  cerr << "                 of width value prior to histogram" << endl;                    
  cerr << "                 estimation." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  irtkGreyImage *ref;
  irtkGreyImage *input;

  int pad = MIN_GREY;
  double sigma = 0;

  irtkGreyPixel *pIn, *pRef;

  int voxels, i, unpaddedCountIn = 0, unpaddedCountRef = 0;
  bool ok;
  float *intsRef;
  float *intsIn;
  float *offsetsIn;
  int offset;

  double temp, scale;

  ref_name    = argv[1];
  argc--;
  argv++;
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  ref = new irtkGreyImage(ref_name);
  input = new irtkGreyImage(input_name);

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
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (sigma > 0){
    // Blur image.
    irtkGaussianBlurringWithPadding<irtkGreyPixel> gaussianBlurring(sigma, pad);
    gaussianBlurring.SetInput (ref);
    gaussianBlurring.SetOutput(ref);
    gaussianBlurring.Run();
    gaussianBlurring.SetInput (input);
    gaussianBlurring.SetOutput(input);
    gaussianBlurring.Run();
  }

  // Count valid voxels.
  voxels = input->GetNumberOfVoxels();
  unpaddedCountIn = 0;
  pIn   = input->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pIn > pad){
      ++unpaddedCountIn;
    }
    ++pIn;
  }

  // Put all valid voxel intensities into an array.
  // Numerical recipes arrays are 1-indexed, so need an extra element at
  // the end.
  intsIn  = new float[unpaddedCountIn];
  offsetsIn = new float[unpaddedCountIn];

  pIn   = input->GetPointerToVoxels();
  unpaddedCountIn = 0;

  for (offset = 0; offset < voxels; ++offset){
    if (*pIn > pad){
      intsIn[unpaddedCountIn] = *pIn;
      offsetsIn[unpaddedCountIn] = offset;
      ++unpaddedCountIn;
    }
    ++pIn;
  }

  sort2(unpaddedCountIn, intsIn - 1, offsetsIn - 1);

  ////////////////

  // Count valid voxels.
  voxels = ref->GetNumberOfVoxels();
  unpaddedCountRef = 0;
  pRef  = ref->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pRef > pad){
      ++unpaddedCountRef;
    }
    ++pRef;
  }

  intsRef  = new float[unpaddedCountRef];

  pRef  = ref->GetPointerToVoxels();
  unpaddedCountRef = 0;

  for (i = 0; i < voxels; ++i){
    if (*pRef > pad){
      intsRef[unpaddedCountRef] = *pRef;
      ++unpaddedCountRef;
    }
    ++pRef;
  }

  sort(unpaddedCountRef, intsRef - 1);

  ////////////////////////////////

  pIn   = input->GetPointerToVoxels();

  scale = (unpaddedCountRef - 1) / ((double) (unpaddedCountIn - 1));

  for (i = 0; i < unpaddedCountIn; ++i){
    temp = floor (i * scale);
    offset = (int) offsetsIn[i];
    *(pIn + offset) = (int) round(intsRef[(int) temp]);
  }

  input->Write(output_name);
}

