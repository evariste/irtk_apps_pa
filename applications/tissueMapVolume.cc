///////////////////////////////////////////////////////////////
// What is the volume of a tissue map?

#include <irtkImage.h>
char *image_name = NULL;

void usage(char *exeName){
  cerr << "usage : " << exeName << " [TissueMap] <-options>" << endl;
  cerr << "TissueMap is a tissue probability map." << endl;
  cerr << "The maximum value in the image is assumed to correspond to a voxel "<< endl;
  cerr << "100% filled with the tissue. This can be overriden with an option." << endl;
  cerr << "Voxels without any tissue present should have a zero value." << endl;
  cerr << "Returns the volume including partially filled voxels." << endl;
  cerr << "Options:" << endl;
  cerr << "-max [value] : Maximum value for full tissue occupancy." << endl;
  cerr << "-t [value]   : Threshold value between 0 (default) and 1. Ignore any " << endl;
  cerr << "               voxels with a fraction of tissue less than [value]." << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 2){
    usage(argv[0]);
  }

  irtkRealPixel min, max;
  double xsize, ysize, zsize, voxelVolume, totalVol;
  double threshold = 0.0;
  int voxels, i;
  bool ok;
  irtkRealPixel *pPix;

  // Parse arguments.
  image_name = argv[1];
  argv++;
  argc--;
  
  irtkRealImage *image = new irtkRealImage(image_name);
  image->GetMinMax(&min, &max);


  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-t") == 0)){
      argc--;
      argv++;
      threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-max") == 0)){
      argc--;
      argv++;
      max = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  image->GetPixelSize(&xsize, &ysize, &zsize);
  voxelVolume = xsize * ysize * zsize;

  pPix     = image->GetPointerToVoxels();
  voxels   = image->GetNumberOfVoxels();
  totalVol = 0.0;

  for (i = 0; i < voxels; ++i){
    if ((*pPix / max) > threshold){
      totalVol += *pPix;
    }
    ++pPix;
  }

  // Convert mm^3 to cm^3
  //cout << 0.001 * totalVol * voxelVolume / max << endl;
  printf("%.3f\n", 0.001 * totalVol * voxelVolume / max);

}
