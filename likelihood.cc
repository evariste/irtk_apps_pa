///////////////////////////////////////////////////////////////
// Given a mean and standard deviation, generate a tissue probability map
// based on the intensities in an image.

#include <irtkImage.h>

#define ROOT_2_PI 2.50662827463100

char *image_name = NULL;
char *output_name = NULL;

void usage(char *exeName){
  cerr << "usage : " << exeName << " [image] [outputLikelihood] [mu] [sigma] <options>" << endl;
  cerr << "Generate likelihoods for a tissue class given the parameters" << endl;
  cerr << "[mu] and [sigma] and the intensities of the input [image]." << endl;
  cerr << "Options: " << endl;
  cerr << "<-pad padding value>" << endl;
  cerr << "<-scale scale value>, default 256" << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 4){
    usage(argv[0]);
  }

  int ok, voxels, i;
  double mu, sigma;
  double pad = -1.0 * FLT_MAX;
  double z, likelihood, denominator;
  double scale = 256;

  image_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;
  mu = atof(argv[1]);
  argv++;
  argc--;
  sigma = atof(argv[1]);
  argv++;
  argc--;

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
    if ((ok == False) && (strcmp(argv[1], "-scale") == 0)){
      argc--;
      argv++;
      scale = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (sigma > 0){
    denominator = 2.0 * sigma * sigma;
  } else {
    cerr << "sigma must be positive ... " << endl;
    exit(1);
  }

  irtkRealImage *image = new irtkRealImage(image_name);
  irtkRealPixel *pPix  = image->GetPointerToVoxels();
  voxels = image->GetNumberOfVoxels();

  double max = -1.0 * FLT_MAX;

  for (i = 0; i < voxels; ++i){
    if (*pPix > pad){
      z = *pPix - mu;
      likelihood = exp( -1.0 * z * z / denominator ) * scale;
      *pPix = likelihood;
      if (*pPix > max)
        max = *pPix;
    } else{
      *pPix = 0;
    }
    ++pPix;
  }

  image->Write(output_name);

}
