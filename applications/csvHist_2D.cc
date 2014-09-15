#include <irtkImage.h>
#include <irtkHistogram.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

char *target_name = NULL, *source_name = NULL;
char *output_name = NULL;
char *dof_name = NULL;

void usage(){
  cout << "csvHist_2D [target] [source] <options>" << endl;

  cout << "Write the histogram of an image to a file in csv format." << endl;
  cout << "Options:" << endl;
  cout << " <-f outputFile>" << endl;
  cout << " <-Tp padding>" << endl;
  cout << " <-bins binCount | -width binWidth>" << endl;
  cout << " <-dofin dof>" << endl;

  cout << "If padding is given, all values <= padding are ignored."  << endl;
  cout << "Only one of -bins or -width should be set, not both." << endl;
  exit(1);
}

int main(int argc, char **argv){

  irtkHistogram_2D<int> hist;

  irtkRealImage target, source;
  irtkTransformation *transformation = NULL;
  irtkInterpolateImageFunction *interpolator = NULL;

  irtkRealPixel tmin = 0, tmax = 0;
  irtkRealPixel smin = 0, smax = 0;

  irtkRealPixel tval, sval;

  int i, j, k, nBins = 0;
  int binx, biny;

  double padding = -1 * FLT_MAX;
  int xdim, ydim, zdim;
  double x, y, z;
  double x1, y1, z1, x2, y2, z2;

  bool ok;
  double width = -1.0;
  double widthX, widthY;

  // Read arguments
  if (argc < 3){
    usage();
  }

  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;


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
    if ((ok == false) && (strcmp(argv[1], "-f") == 0)){
      argc--;
      argv++;
      output_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      padding = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bins") == 0)){
      argc--;
      argv++;
      nBins = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-width") == 0)){
      argc--;
      argv++;
      width = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    transformation = new irtkRigidTransformation;
  }


  // Initialise histogram.
  target.Read(target_name);
  source.Read(source_name);


  interpolator = new irtkLinearInterpolateImageFunction;
  interpolator->SetInput(&source);
  interpolator->Initialize();
  // Calculate the source image domain in which we can interpolate
  interpolator->Inside(x1, y1, z1, x2, y2, z2);

  xdim = target.GetX();
  ydim = target.GetY();
  zdim = target.GetZ();

  tmin = FLT_MAX;
  smin = FLT_MAX;

  tmax = -1.0 * FLT_MAX;
  smax = -1.0 * FLT_MAX;

  // What are the min and max values corresponding to the unpadded region
  // in the target?
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){

        x = i;
        y = j;
        z = k;

        tval = target.Get(i, j, k);

        if (tval > padding){

          target.ImageToWorld(x, y, z);
          transformation->Transform(x, y, z);
          source.WorldToImage(x, y, z);

          if ((x > x1) && (x < x2) &&
              (y > y1) && (y < y2) &&
              (z > z1) && (z < z2)){

            if (tval >  tmax)
              tmax = tval;
            if (tval < tmin)
              tmin = tval;

            sval = interpolator->EvaluateInside(x, y, z);

            if (sval >  smax)
              smax = sval;
            if (sval < smin)
              smin = sval;
          }

        }
      }
    }
  }

  // Set up the histogram.

  cerr << "tmax " << tmax << " tmin " << tmin << endl;
  cerr << "smax " << smax << " smin " << smin << endl;

  hist.PutMin(floor(tmin), floor(smin));
  hist.PutMax(floor(tmax) + 1, floor(smax) + 1);

  if (width > 0){
    if (nBins > 0){
      // Both width and bins have been set.
      cerr << "Only one of -bins or -width should be set" << endl;
      exit(1);
    } else {
      // Width only set.
      hist.PutWidth(width, width);
      hist.GetNumberOfBins(&binx, &biny);
      cerr << "Histogram assigned width " << width << ", using " << binx;
      cerr << " x " << biny << " bins." << endl;
    }
  } else {
    // Width not set.

    if (nBins > 0){
      // Bins only set.
      hist.PutNumberOfBins(nBins, nBins);
      hist.GetWidth(&widthX, &widthY);
      cerr << "Histogram assigned " << nBins << " bins, using widths ";
      cerr << widthX << " x " << widthY << endl;
    } else {
      // Neither set.
      hist.GetWidth(&widthX, &widthY);
      hist.GetNumberOfBins(&binx, &biny);
      cerr << "Histogram uses " << binx << " x " << biny  << " bins with widths ";
      cerr << widthX << " x " << widthY << endl;
    }
  }

  // Populate histogram.
  hist.Reset();
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        x = i;
        y = j;
        z = k;

        tval = target.Get(i, j, k);
        if (tval > padding){

          target.ImageToWorld(x, y, z);
          transformation->Transform(x, y, z);
          source.WorldToImage(x, y, z);

          if ((x > x1) && (x < x2) &&
              (y > y1) && (y < y2) &&
              (z > z1) && (z < z2)){

            sval = interpolator->EvaluateInside(x, y, z);

            hist.AddSample(tval, sval);
          }
        }
      }
    }
  }


  cerr << hist.NumberOfSamples() << endl;


  if (output_name != NULL){
    // Write histogram.
    ofstream fileOut(output_name);

    if(!fileOut){
      cerr << "Can't open file " << output_name << endl;
      exit(1);
    }

    hist.GetNumberOfBins(&binx, &biny);

    for (i = 0; i < binx - 1; ++i){
      fileOut << hist.BinToValX(i) << ",";
    }
    fileOut << hist.BinToValX(i) << endl;


    for (j = 0; j < biny - 1; ++j){
      fileOut << hist.BinToValY(j) << ",";
    }
    fileOut << hist.BinToValY(j) << endl;


    for (j = 0; j < biny; ++j){
      for (i = 0; i < binx - 1; ++i){
        fileOut << hist(i, j) << ",";
      }
      // Avoid trailing comma.
      fileOut << hist(i, j) << endl;
    }
    fileOut.close();

  } else {

    // write to std out.
    for (i = 0; i < binx - 1; ++i){
      cout << hist.BinToValX(i) << ",";
    }
    cout << hist.BinToValX(i) << endl;

    for (j = 0; j < biny - 1; ++j){
      cout << hist.BinToValY(j) << ",";
    }
    cout << hist.BinToValY(j) << endl;

    for (j = 0; j < biny; ++j){
      for (i = 0; i < binx - 1; ++i){
        cout << hist(i, j) << ",";
      }
      // Avoid trailing comma.
      cout << hist(i, j) << endl;
    }
  }


}
