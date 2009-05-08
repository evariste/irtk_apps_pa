///////////////////////////////////////////////////////////////
// What is the volume of a tissue map?

#include <irtkImage.h>
char *image_name = NULL;

void usage(char *exeName){
  cerr << "usage : " << exeName << " [TissueMap] [maxValue] <-options>" << endl;
  cerr << "TissueMap is a tissue probability map." << endl;
  cerr << "maxValue is the value given when a voxel is 100% filled with the tissue." << endl;
  cerr << "Voxels without any tissue present should have a zero value." << endl;
  cerr << "Returns the volume including partially filled voxels." << endl;
  cerr << "Options:" << endl;
  cerr << "-t [value] : Threshold value between 0 and 1. Ignore any voxels" << endl;
  cerr << "             with a fraction of tissue less than [value]." << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 3){
    usage(argv[0]);
  }

  irtkRealPixel max;
  double xsize, ysize, zsize, voxelVolume, totalVol;
  double threshold = -1.0 * FLT_MAX;
  int voxels, i, ok;
  irtkRealPixel *pPix;

  // Parse arguments.
  image_name = argv[1];
  argv++;
  argc--;
  max = atof(argv[1]);
  argv++;
  argc--;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-t") == 0)){
      argc--;
      argv++;
      threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  /////////////////////////////////////////

  irtkRealImage *image = new irtkRealImage(image_name);

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
  cout << 0.001 * totalVol * voxelVolume / max << endl;
}
