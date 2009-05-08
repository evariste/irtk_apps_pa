#include <irtkImage.h>


int main(int argc, char* argv[])
{
  int intersectionCount, count1, count2, unionCount;
  double overlap, si;
  double threshold = 500;
  int ok;
  int reportOverlap = True;

  if (argc < 3){
    cerr << "Usage: overlap [image1] [image2] <-threshold> <-si>" << endl;
    cerr << "-threshold = value above which a voxel is considered as labeled by either image." << endl;
    cerr << "Default threshold = 500" << endl;
    cerr << "If SI is required instead, use -si flag." << endl;
    exit(0);
  }

  // Read first image
  irtkGreyImage image1(argv[1]);
  argc--;      argv++;
  // Read second image  
  irtkGreyImage image2(argv[1]);
  argc--;      argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;      argv++;
      threshold = atof(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-si") == 0)){
      reportOverlap = False;
      argc--;      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (image1.GetX() != image2.GetX() || image1.GetY() != image2.GetY() || image1.GetZ() != image2.GetZ()){
    cerr << "Images must have the same numbers of voxels along each dimension." << endl;
    exit(1);
  }

  double x1, y1, z1, x2, y2, z2;
  double errx, erry, errz;

  image1.GetPixelSize(&x1, &y1, &z1);
  image2.GetPixelSize(&x2, &y2, &z2);

  errx = fabs(x1 - x2) / x1;
  erry = fabs(y1 - y2) / y1;
  errz = fabs(z1 - z2) / z1;


  if (errz > 0.01 || erry > 0.01 || errz > 0.01){
    cerr << "overlap.cc : Warning : differences between voxel sizes detected." << endl;
  }


  irtkGreyPixel *ptrPix1 = image1.GetPointerToVoxels();
  irtkGreyPixel *ptrPix2 = image2.GetPointerToVoxels();

  count1 = 0;
  count2 = 0;
  intersectionCount = 0;

  for (int i=0; i< image1.GetNumberOfVoxels(); ++i){
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
    cout << overlap << "";
  } else {
    cout << si << "";
  }
}
