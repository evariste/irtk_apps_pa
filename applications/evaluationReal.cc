#include <irtkImage.h>

#include <irtkImageFunction.h>

#include <irtkHistogram.h>

#include <irtkTransformation.h>

// Default filenames
char *source_name = NULL, *target_name, *dof_name  = NULL;

// Default number of bins
int nbins_x = 0, nbins_y = 0;

void usage()
{
  cerr << "Usage: evaluationReal [target] [source] <options>\n" << endl;
  cerr << "where <options> is one or more of the following: \n" << endl;
  cerr << "<-dofin file>      Transformation" << endl;
  cerr << "<-nbins_x no>      Number of bins for target" << endl;
  cerr << "<-nbins_y no>      Number of bins for source" << endl;
  cerr << "<-Rx1 pixel>       Region of interest" << endl;
  cerr << "<-Ry1 pixel>       Region of interest" << endl;
  cerr << "<-Rz1 pixel>       Region of interest" << endl;
  cerr << "<-Rx2 pixel>       Region of interest" << endl;
  cerr << "<-Ry2 pixel>       Region of interest" << endl;
  cerr << "<-Rz2 pixel>       Region of interest" << endl;
  cerr << "<-Tp  value>       Padding value in target" << endl;
  cerr << "<-linear>          Linear interpolation" << endl;
  cerr << "<-bspline>         B-spline interpolation" << endl;
  cerr << "<-cspline>         Cubic spline interpolation" << endl;
  cerr << "<-sinc>            Sinc interpolation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation *transformation = NULL;
  irtkInterpolateImageFunction *interpolator = NULL;
  irtkRealPixel target_min, source_min, target_max, source_max;
  bool ok;
  int x, y, z, t, i1, j1, k1, i2, j2, k2;
  double x1, y1, z1, x2, y2, z2, Tp, widthx, widthy, val;

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

  // Fix ROI
  Tp = -1.0 * FLT_MAX;
  i1 = 0;
  j1 = 0;
  k1 = 0;
  i2 = target.GetX();
  j2 = target.GetY();
  k2 = target.GetZ();

  // Fix no. of bins;
  nbins_x = 0;
  nbins_y = 0;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nbins_x") == 0)){
      argc--;
      argv++;
      nbins_x = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nbins_y") == 0)){
      argc--;
      argv++;
      nbins_y = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      Tp = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)){
      argc--;
      argv++;
      i1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)){
      argc--;
      argv++;
      i2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)){
      argc--;
      argv++;
      j1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)){
      argc--;
      argv++;
      j2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)){
      argc--;
      argv++;
      k1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)){
      argc--;
      argv++;
      k2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)){
      argc--;
      argv++;
      interpolator = new irtkLinearInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)){
      argc--;
      argv++;
      interpolator = new irtkBSplineInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)){
      argc--;
      argv++;
      interpolator = new irtkCSplineInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)){
      argc--;
      argv++;
      interpolator = new irtkSincInterpolateImageFunction;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
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
  cout << "Min and max of X is " << target_min 
       << " and " << target_max << endl;
  cout << "Min and max of Y is " << source_min 
       << " and " << source_max << endl;

  // Calculate number of bins to use
  if (nbins_x == 0){
    nbins_x = (int) round(target_max - target_min) + 1;
  }

  if (nbins_y == 0){
    nbins_y = (int) round(source_max - source_min) + 1;
  }

  // Create default interpolator if necessary
  if (interpolator == NULL){
    interpolator = 
      new irtkNearestNeighborInterpolateImageFunction;
  }
  interpolator->SetInput(&source);
  interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  interpolator->Inside(x1, y1, z1, x2, y2, z2);

  // Create histogram
  irtkHistogram_2D<int> histogram(nbins_x, nbins_y);

  widthx = (target_max - target_min) / (nbins_x - 1.0);
  widthy = (source_max - source_min) / (nbins_y - 1.0);

  histogram.PutMin(target_min - 0.5*widthx, source_min - 0.5*widthy);
  histogram.PutMax(target_max + 0.5*widthx, source_max + 0.5*widthy);

  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    transformation = new irtkRigidTransformation;
  }

  target_min = FLT_MAX;
  source_min = FLT_MAX;
  target_max = -1.0 * FLT_MAX;
  source_max = -1.0 * FLT_MAX;

   
  for (t = 0; t < target.GetT(); t++) {

    cout << "Timepoint " << t+1 << endl;

    histogram.Reset();

    // Fill histogram
  for (z = 0; z < target.GetZ(); z++){
    for (y = 0; y < target.GetY(); y++){
      for (x = 0; x < target.GetX(); x++){

        val = target(x, y, z, t);

	if (val > Tp){

          if (val > target_max)
            target_max = val;
          if (val < target_min)
            target_min = val;

	  irtkPoint p(x, y, z);
	  // Transform point into world coordinates
	  target.ImageToWorld(p);
	  // Transform point
	  transformation->Transform(p);
	  // Transform point into image coordinates
	  source.WorldToImage(p);
	  if ((p._x > x1) && (p._x < x2) &&
	      (p._y > y1) && (p._y < y2) &&
	      (p._z > z1) && (p._z < z2)){

            val = interpolator->EvaluateInside(p._x, p._y, p._z, t);

            histogram.AddSample(target(x, y, z, t), val);
            if (val >  source_max)
              source_max = val;
            if (val < source_min)
              source_min = val;
	  }
	} 
      }
    }
  }

  cout << "ROI Min and max of X is " << target_min 
       << " and " << target_max << endl;
  cout << "ROI Min and max of Y is " << source_min 
       << " and " << source_max << endl;
  int binx, biny;
  histogram.GetNumberOfBins(&binx, &biny);
  cout << "Number of bins  X x Y : " << binx << " x " << biny << endl;
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
  }

}

