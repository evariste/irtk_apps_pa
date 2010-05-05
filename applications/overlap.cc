#include <irtkImage.h>


int main(int argc, char* argv[])
{
  int intersectionCount, count1, count2, unionCount;
  double overlap, si;
  double threshold;
  int ok;
  int reportOverlap = True;
  double x1, y1, z1, x2, y2, z2;
  double errx, erry, errz;
  float min1, max1, min2, max2;
  int fuzzy = False;
  double sigmaMin = 0.0;
  double sigmaMax = 0.0;
  double sigmaMinForVol, sigmaMaxForVol;
  double *labelVols;
  int i, t, noOfVoxels, voxelsPerVolume;

  if (argc < 3){
    cerr << "Usage: overlap [image1] [image2] <-threshold> <-si>" << endl;
    cerr << "Provide overlap of a label in two images (Tanimoto)" << endl;

    cerr << "-threshold   Value above which a voxel is considered as labeled in either image." << endl;
    cerr << "             Default threshold = half maximum value across both images." << endl;
    cerr << "-si          Give similarity index (Dice coefficient) instead of Tanimoto overlap." << endl;
    cerr << "-fuzzy       Give fuzzy overlap, the voxelwise sum of the minimum values over the " << endl;
    cerr << "             the voxelwise sum of maximum values.  Previous options no longer apply." << endl;
    exit(0);
  }

  // Read first image
  irtkRealImage image1(argv[1]);
  argc--;
  argv++;

  // Read second image
  irtkRealImage image2(argv[1]);
  argc--;
  argv++;

  image1.GetMinMax(&min1, &max1);
  image2.GetMinMax(&min2, &max2);
  max1 = (max1 > max2) ? max1 : max2;
  threshold = 0.5 * max1;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;
      threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-si") == 0)){
      reportOverlap = False;
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-fuzzy") == 0)){
      fuzzy = True;
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (image1.GetX() != image2.GetX() ||
      image1.GetY() != image2.GetY() ||
      image1.GetZ() != image2.GetZ() ||
      image1.GetT() != image2.GetT()){
    cerr << "Images must have the same dimensions." << endl;
    exit(1);
  }


  image1.GetPixelSize(&x1, &y1, &z1);
  image2.GetPixelSize(&x2, &y2, &z2);

  errx = fabs(x1 - x2) / x1;
  erry = fabs(y1 - y2) / y1;
  errz = fabs(z1 - z2) / z1;


  if (errz > 0.01 || erry > 0.01 || errz > 0.01){
    cerr << "overlap : Warning : differences between voxel sizes detected." << endl;
  }



  irtkRealPixel *ptrPix1, *ptrPix2;

  noOfVoxels = image1.GetNumberOfVoxels();


  if (fuzzy == True){

    double *labelVols = new double[image1.GetT()];
    voxelsPerVolume = noOfVoxels / image1.GetT();

    for (t = 0; t < image1.GetT(); t++){
      ptrPix1 = image1.GetPointerToVoxels(0, 0, 0, t);
      ptrPix2 = image2.GetPointerToVoxels(0, 0, 0, t);

      for (i = 0; i < voxelsPerVolume; ++i){
        labelVols[t] += *ptrPix1;
        labelVols[t] += *ptrPix2;
        ++ptrPix1;
        ++ptrPix2;
      }
    }

    for (t = 0; t < image1.GetT(); t++){
      if (labelVols[t] <= 0){
        cerr << "Zero volume for label (index t = " << t << ")" << endl;
        exit(1);
      }
    }

    for (t = 0; t < image1.GetT(); t++){
      ptrPix1 = image1.GetPointerToVoxels(0, 0, 0, t);
      ptrPix2 = image2.GetPointerToVoxels(0, 0, 0, t);

      sigmaMinForVol = 0.0;
      sigmaMaxForVol = 0.0;
      for (i = 0; i < voxelsPerVolume; ++i){
        sigmaMinForVol += (*ptrPix1 < *ptrPix2) ? *ptrPix1 : *ptrPix2;
        sigmaMaxForVol += (*ptrPix1 > *ptrPix2) ? *ptrPix1 : *ptrPix2;
        ++ptrPix1;
        ++ptrPix2;
      }

      // Weight inversely by volume (Crum TMI 2006)
      sigmaMin += sigmaMinForVol / labelVols[t];
      sigmaMax += sigmaMaxForVol / labelVols[t];
    }


    cout << sigmaMin / sigmaMax << endl;

  }else{

    // Plain hard labels.

    ptrPix1 = image1.GetPointerToVoxels();
    ptrPix2 = image2.GetPointerToVoxels();

    count1 = 0;
    count2 = 0;
    intersectionCount = 0;

    for (i=0; i< noOfVoxels; ++i){
      if ((*ptrPix1 < threshold) && (*ptrPix2 < threshold)){
        ++ptrPix1;
        ++ptrPix2;
        continue;
      }

      if ((*ptrPix1 >= threshold) && (*ptrPix2 >= threshold)){
        intersectionCount++;
      }

      if (*ptrPix1 >= threshold){
        count1++;
      }

      if (*ptrPix2 >= threshold){
        count2++;
      }

      ++ptrPix1;
      ++ptrPix2;
    }

    unionCount = count1 + count2 - intersectionCount;
    overlap = (double) intersectionCount / (double)unionCount;
    si = 2.0 * (double) intersectionCount / (double) (count1 + count2);

    if (reportOverlap){
      cout << overlap << endl;
    } else {
      cout << si << endl;
    }

  }


}
