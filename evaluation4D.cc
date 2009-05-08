/*=========================================================================
 
  Library   : packages/applications
  Module    : $RCSfile: evaluation.cc,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2000-2001
  Purpose   :
  Date      : $Date: 2005-12-22 15:20:06 $
  Version   : $Revision: 1.7 $
  Changes   : $Locker:  $
	      $Log: evaluation.cc,v $
	      Revision 1.7  2005-12-22 15:20:06  dr
	      Updated to reflect changes to transformation classes

	      Revision 1.6  2005/02/01 10:29:55  dr
	      Updated copyright notices
	
	      Revision 1.5  2003/10/28 18:52:18  dr
	      Nearest neighbor interpolation is default
	
	      Revision 1.4  2003/10/07 10:12:42  dr
	      Image are now always interpolated
	
	      Different interpolation schemes have been implemented
	
	      Revision 1.3  2003/08/18 13:02:39  dr
	      Modified interpolation framework
	
	      Revision 1.2  2003/04/17 14:10:51  dr
	      Merged branch
	
	      Revision 1.1.1.1.2.2  2002/09/15 11:36:31  dr
	      Label consistency and Kappa statistic are only calculated
	      if number of histogram bins is identical
	
	      Revision 1.1.1.1.2.1  2001/11/26 15:55:24  dr
	      Added changes from Alex
	
	      Revision 1.1.1.1  2000/10/29 15:12:16  dr
	      Imported sources
	

=========================================================================*/

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
  cerr << "Usage: evaluation [target] [source] <options>\n" << endl;
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
  irtkGreyPixel target_min, source_min, target_max, source_max;
  int ok, x, y, z, i1, j1, k1, i2, j2, k2, Tp; 
  double x1, y1, z1, x2, y2, z2; 
  int t;

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
  irtkGreyImage target(target_name);
  cout << "done" << endl;
  // Read source image
  cout << "Reading source ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

  // Fix ROI
  Tp = -INT_MAX;
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
      Tp = atoi(argv[1]);
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
    if ((ok == False) && (strcmp(argv[1], "-linear") == 0)){
      argc--;
      argv++;
      interpolator = new irtkLinearInterpolateImageFunction;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bspline") == 0)){
      argc--;
      argv++;
      interpolator = new irtkBSplineInterpolateImageFunction;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-cspline") == 0)){
      argc--;
      argv++;
      interpolator = new irtkCSplineInterpolateImageFunction;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sinc") == 0)){
      argc--;
      argv++;
      interpolator = new irtkSincInterpolateImageFunction;
      ok = True;
    }
    if (ok == False){
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
  if (nbins_x == 0) nbins_x = target_max - target_min + 1;
  if (nbins_y == 0) nbins_y = source_max - source_min + 1;

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
  irtkHistogram_2D histogram(nbins_x, nbins_y);
  histogram.PutMin(target_min - 0.5, source_min - 0.5);
  histogram.PutMax(target_max + 0.5, source_max + 0.5);

  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    transformation = new irtkRigidTransformation;
  }
   

  // Check t dimensions equal.
  if (target.GetT() != source.GetT()){
    cerr << "Images have different T dimensions." << endl;
    exit(1);
  }

  for (t = 0; t < target.GetT(); t++){

    cout << "Timepoint " << t+1 << endl;

    histogram.Reset();

    // Fill histogram
    for (z = 0; z < target.GetZ(); z++){
      for (y = 0; y < target.GetY(); y++){
        for (x = 0; x < target.GetX(); x++){
          if (target(x,y,z,t) > Tp){
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
      histogram.AddSample(target(x, y, z, t), round(interpolator->EvaluateInside(p._x, p._y, p._z, t)));
            }
          } 
        }
      }
    }


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

    cout << endl;

  }

  cout << "------------------------" << endl;
  cout << "All timepoints combined:" << endl;
  cout << "------------------------" << endl;
  
  histogram.Reset();

  for (t = 0; t < target.GetT(); t++){
    for (z = 0; z < target.GetZ(); z++){
      for (y = 0; y < target.GetY(); y++){
        for (x = 0; x < target.GetX(); x++){

          if (target(x,y,z,t) > Tp){
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
              histogram.AddSample(target(x, y, z, t),
                                  round(interpolator->EvaluateInside(p._x, p._y, p._z, t)));
            }
          } 
        }
      }
    }
  }

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
  
  cout << endl;

}

