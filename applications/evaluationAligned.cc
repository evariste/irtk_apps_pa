#include <irtkImage.h>

#include <irtkImageFunction.h>

#include <irtkHistogram.h>

#include <irtkTransformation.h>

// Default filenames
char *source_name = NULL, *target_name, *dof_name  = NULL;
char *mask_name = NULL;

// Default number of bins
int nbins_x = 256, nbins_y = 256;

void usage()
{
  cerr << "Usage: evaluationAligned [target] [source] <options>\n" << endl;
  cerr << "[target] [source] are assumed to be aligned and on the same lattice." << endl << endl;

  cerr << "where <options> is one or more of the following: \n" << endl;
  cerr << "<-mask file>       Mask for region to be evaluated." << endl;
  cerr << "<-nbins_x no>      Number of bins for target (default 256)" << endl;
  cerr << "<-nbins_y no>      Number of bins for source (default 256)" << endl;
  cerr << "<-Rx1 pixel>       Region of interest" << endl;
  cerr << "<-Ry1 pixel>       Region of interest" << endl;
  cerr << "<-Rz1 pixel>       Region of interest" << endl;
  cerr << "<-Rx2 pixel>       Region of interest" << endl;
  cerr << "<-Ry2 pixel>       Region of interest" << endl;
  cerr << "<-Rz2 pixel>       Region of interest" << endl;
  cerr << "<-Tp  value>       Padding value in target" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkRealPixel target_min, source_min, target_max, source_max;
  int ok, x, y, z, i1, j1, k1, i2, j2, k2; 
  double x1, y1, z1, x2, y2, z2, Tp, widthx, widthy, val;
  double unbinnedSSD;
  int n;

  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkRealImage target(target_name);
  cout << "done" << endl;
  // Read source image
  cout << "Reading source ... "; cout.flush();
  irtkRealImage source(source_name);
  cout << "done" << endl;

  if (target.GetX() != source.GetX() || target.GetY() != source.GetY() || target.GetZ() != source.GetZ()){
    cerr << "Images are not on the same lattice\n" << endl;
    exit(1);
  }

  // Fix ROI
  Tp = -1.0 * FLT_MAX;
  i1 = 0;
  j1 = 0;
  k1 = 0;
  i2 = target.GetX();
  j2 = target.GetY();
  k2 = target.GetZ();

  // Fix no. of bins;
  nbins_x = 256;
  nbins_y = 256;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-nbins_x") == 0)){
      argc--;
      argv++;
      nbins_x = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-nbins_y") == 0)){
      argc--;
      argv++;
      nbins_y = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      Tp = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)){
      argc--;
      argv++;
      i1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)){
      argc--;
      argv++;
      i2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)){
      argc--;
      argv++;
      j1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)){
      argc--;
      argv++;
      j2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)){
      argc--;
      argv++;
      k1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)){
      argc--;
      argv++;
      k2 = atoi(argv[1]);
      argc--;
      argv++;
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
      usage();
    }
  } 


  irtkGreyImage mask;
  // Is there a mask to use?
  if (mask_name != NULL){
    mask.Read(mask_name);
    if (target.GetX() != mask.GetX() || target.GetY() != mask.GetY() || target.GetZ() != mask.GetZ()){
      cerr << "Mask is not on the same lattice\n" << endl;
      exit(1);
    }
  } else {
    mask.Read(target_name);
    irtkGreyPixel *ptr2mask = mask.GetPointerToVoxels();
    for (int i = 0; i < mask.GetNumberOfVoxels(); ++i){
      *ptr2mask = 1;
      ++ptr2mask;
    }
  }

  // If there is an region of interest, use it
  if ((i1 != 0) || (i2 != target.GetX()) ||
      (j1 != 0) || (j2 != target.GetY()) || 
      (k1 != 0) || (k2 != target.GetZ())){
    target = target.GetRegion(i1, j1, k1, i2, j2, k2);
    source = source.GetRegion(i1, j1, k1, i2, j2, k2);
  }

  // Set min and max of histogram
  target.GetMinMax(&target_min, &target_max);
  source.GetMinMax(&source_min, &source_max);

  // If user has set bins to zero, assume width 1 bins required.
  if (nbins_x == 0){
    nbins_x = (int) round(target_max - target_min) + 1;
  }

  if (nbins_y == 0){
    nbins_y = (int) round(source_max - source_min) + 1;
  }


  cout << "Min and max of X is " << target_min 
       << " and " << target_max << endl;
  cout << "Min and max of Y is " << source_min 
       << " and " << source_max << endl;


  // Create histogram
  irtkHistogram_2D histogram(nbins_x, nbins_y);

  widthx = (target_max - target_min) / (nbins_x - 1.0);
  widthy = (source_max - source_min) / (nbins_y - 1.0);

  histogram.PutMin(target_min - 0.5*widthx, source_min - 0.5*widthy);
  histogram.PutMax(target_max + 0.5*widthx, source_max + 0.5*widthy);


  target_min = FLT_MAX;
  source_min = FLT_MAX;
  target_max = -1.0 * FLT_MAX;
  source_max = -1.0 * FLT_MAX;
  int maskVal;

  // Fill histogram
  for (z = 0; z < target.GetZ(); z++){
    for (y = 0; y < target.GetY(); y++){
      for (x = 0; x < target.GetX(); x++){

        val = target(x, y, z);
        maskVal = mask(x, y, z);

	if (val > Tp && maskVal > 0){

          if (val > target_max)
            target_max = val;
          if (val < target_min)
            target_min = val;

          histogram.AddSample(target(x, y, z), source(x, y, z));

          if (source(x, y, z) >  source_max)
            source_max = source(x, y, z);
          if (source(x, y, z) < source_min)
            source_min = source(x, y, z);

          unbinnedSSD += (target(x, y, z) - source(x, y, z)) * (target(x, y, z) - source(x, y, z));
          ++n; 

	} 
      }
    }
  }

  cout << "ROI Min and max of X is " << target_min 
       << " and " << target_max << endl;
  cout << "ROI Min and max of Y is " << source_min 
       << " and " << source_max << endl;
  cout << "Number of bins  X x Y : " << histogram.GetNumberOfBinsX() << " x " << histogram.GetNumberOfBinsY() << endl;
  cout << "Number of Samples: "     << histogram.NumberOfSamples() << endl;
  cout << "Mean of X: "             << histogram.MeanX() << endl;
  cout << "Mean of Y: "             << histogram.MeanY() << endl;
  cout << "Variance of X: "         << histogram.VarianceX() << endl;
  cout << "Variance of Y: "         << histogram.VarianceY() << endl;
  cout << "Covariance: "            << histogram.Covariance() << endl;
  cout << "Crosscorrelation: "      << histogram.CrossCorrelation() << endl;
  cout << "Sums of squared diff.: " << histogram.SumsOfSquaredDifferences() / (double)histogram.NumberOfSamples() << endl;
  cout << "Marginal entropy of X: " << histogram.EntropyX() << endl;
  cout << "Marginal entropy of Y: " << histogram.EntropyY() << endl;
  cout << "Joint entropy: "         << histogram.JointEntropy() << endl;
  cout << "Mutual Information: "    << histogram.MutualInformation() << endl;
  cout << "Norm. Mutual Information: " << histogram.NormalizedMutualInformation() << endl;
  cout << "Correlation ratio C(X|Y): " << histogram.CorrelationRatioXY() << endl;
  cout << "Correlation ratio C(Y|X): " << histogram.CorrelationRatioYX() << endl;
  if (nbins_x == nbins_y){
    cout << "Label consistency: " << histogram.LabelConsistency() << endl;
    cout << "Kappa statistic: " << histogram.Kappa() << endl;
  }
  cout << "Unbinned SSD :            " << unbinnedSSD / n << endl;

}

