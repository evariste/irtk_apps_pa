///////////////////////////////////////////////////////////////
// Stub to test the irtkParzenHistogram_2D histogram class

#include <sys/types.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <irtkImage.h>
#include <nr.h>

#include <irtkHistogram.h>

#include <irtkImageFunction.h>

#include <irtkTransformation.h>

char *target_name = NULL, *source_name = NULL, *dof_name  = NULL;

void usage(){
  cout << "parzenStub.cc [target] [source] <options>" << endl;
  cout << "Options:" << endl;
  cout << "-dofin       [filename]" << endl;
  cout << "-Tp          [padding value]" << endl;
  cout << "-kernel      [kernel width]" << endl;
  cout << "-kernelX     [kernel width X]" << endl;
  cout << "-kernelY     [kernel width Y]" << endl;
  cout << "-bins        [no.]" << endl;
  cout << "-integralPts [no.]" << endl;
  cout << "" << endl;
  exit(0);
}

int main(int argc, char **argv){

  int ok, i, j, k, xdim, ydim, zdim;
  int integralPoints = 100, bins = 100;

  double x, y, z, x1, x2, y1, y2, z1, z2;
  double kernelWidthX = 10, kernelWidthY = 10;

  irtkRealPixel tgtMin, tgtMax, srcMin, srcMax;
  irtkRealPixel padding, targetVal, sourceVal;

  if (argc < 3){
    usage();
  }

  // Stuff for rand calls.
  //timeval tv;
  //gettimeofday(&tv, NULL);
  //long init = tv.tv_usec;
  //long temp = -1*init;
  //(void) ran2(&temp);

  time_t seconds;
  long ran2Seed, ran2initialSeed;

  seconds = time(NULL);
  ran2Seed = seconds;
  ran2initialSeed = -1 * ran2Seed;
  (void) ran2(&ran2initialSeed);

  // Read input.
  target_name = argv[1];
  argv++;
  argc--;
  source_name = argv[1];
  argv++;
  argc--;


  irtkRealImage *target = new irtkRealImage(target_name);
  irtkRealImage *source = new irtkRealImage(source_name);

  irtkTransformation *transformation = NULL;

  irtkInterpolateImageFunction *interpolator
    = new irtkLinearInterpolateImageFunction;

  xdim = target->GetX();
  ydim = target->GetY();
  zdim = target->GetZ();

  padding = -1.0 * FLT_MAX;


  // Parse args.
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
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      padding = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-kernel") == 0)){
      argc--;
      argv++;
      kernelWidthX = kernelWidthY = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-kernelX") == 0)){
      argc--;
      argv++;
      kernelWidthX = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-kernelY") == 0)){
      argc--;
      argv++;
      kernelWidthY = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bins") == 0)){
      argc--;
      argv++;
      bins = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-integralPts") == 0)){
      argc--;
      argv++;
      integralPoints = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    transformation = new irtkRigidTransformation;
  }

  interpolator->SetInput(source);
  interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  interpolator->Inside(x1, y1, z1, x2, y2, z2);

  /////////////////////////////////////////////////////////
  // mins and maxs.

  target->GetMinMax(&tgtMin, &tgtMax);
  source->GetMinMax(&srcMin, &srcMax);

  // Is there a genuine padding value?
  if (padding > -1.0 * FLT_MAX){
    // Find the min max among unpadded voxels.
    tgtMin = FLT_MAX;
    srcMin = FLT_MAX;
    tgtMax = -1.0 * FLT_MAX;
    srcMax = -1.0 * FLT_MAX;

    for (k = 0; k < zdim; ++k){
      for (j = 0; j < ydim; ++j){
        for (i = 0; i < xdim; ++i){
          x = i;
          y = j;
          z = k;
          targetVal = target->Get(i, j, k);

          if (targetVal <= padding)
            continue;

          target->ImageToWorld(x, y, z);
          transformation->Transform(x, y, z);
          source->WorldToImage(x, y, z);

          if ((x > x1) && (x < x2) &&
              (y > y1) && (y < y2) &&
              (z > z1) && (z < z2)){

            if (targetVal < tgtMin)
              tgtMin = targetVal;

            if (targetVal > tgtMax)
              tgtMax = targetVal;

            sourceVal = interpolator->EvaluateInside(x, y, z);

            if (sourceVal < srcMin)
              srcMin = sourceVal;

            if (sourceVal > srcMax)
              srcMax = sourceVal;

          }
        }
      }
    }
  } // if (genuine padding)


  tgtMin = floor(tgtMin);
  srcMin = floor(srcMin);
  tgtMax = 1 + floor(tgtMax);
  srcMax = 1 + floor(srcMax);

  cout << "Using target min and max : \t"  << tgtMin << "\t" << tgtMax << endl;
  cout << "Using source min and max : \t"  << srcMin << "\t" << srcMax << endl;


  /////////////////////////////////////////////////////////
  // The main bit.

  irtkParzenHistogram_2D parzenHist(bins, bins);

  parzenHist.PutMin(floor(tgtMin)    , floor(srcMin));
  parzenHist.PutMax(1 + floor(tgtMax), 1 + floor(srcMax));

  parzenHist.PutIntegralPoints(integralPoints);
  parzenHist.PutKernelWidthX(kernelWidthX);
  parzenHist.PutKernelWidthY(kernelWidthY);
  parzenHist.Initialize();

  // Fill the underlying histogram.
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        x = i;
        y = j;
        z = k;
        targetVal = target->Get(i, j, k);

        if (targetVal <= padding)
          continue;

        target->ImageToWorld(x, y, z);
        transformation->Transform(x, y, z);
        source->WorldToImage(x, y, z);

        if ((x > x1) && (x < x2) &&
            (y > y1) && (y < y2) &&
            (z > z1) && (z < z2)){

            parzenHist.AddSample(targetVal, interpolator->EvaluateInside(x, y, z));

        }
      }
    }
  }

  ////////////////////////////////////////////////
  // Bunch of reporting:
  cout << "Integral        " << parzenHist.Integrate() << endl;

  parzenHist.Print("parzen.hist");
  parzenHist.PrintMarginals("parzen.marginals");

  parzenHist.PrintKernel("parzen.kernel");
  parzenHist.PrintDerivativeX("parzen.dx");
  parzenHist.PrintDerivativeY("parzen.dy");

  // Write the histogram to an image.
  double val;
  double xaxis[3], yaxis[3], zaxis[3];
  int srcDim, tgtDim;
  double aspect, dTgt, dSrc;
  irtkPoint origin;

  tgtDim = 100;
  aspect = (tgtMax - tgtMin) / ((double) (srcMax - srcMin));
  srcDim = (int) (round(tgtDim * aspect));

  dTgt = ((double) tgtMax - tgtMin) / ((double) tgtDim - 1);
  dSrc = ((double) srcMax - srcMin) / ((double) srcDim - 1);

  origin._x = ((double) tgtMin + tgtMax ) / 2.0;
  origin._y = ((double) srcMin + srcMax ) / 2.0;
  origin._z = 0;

  for (i = 0; i < 3; ++i){
    xaxis[i] = yaxis[i] = zaxis[i] = 0;
  }
  xaxis[0] = yaxis[1] = zaxis[2] = 1.0;

  irtkImageAttributes attr;
  attr._x = tgtDim;
  attr._y = srcDim;
  attr._z = 1;
  attr._dx = dTgt;
  attr._dy = dSrc;
  attr._dz = 1;
  attr._xorigin = origin._x;
  attr._yorigin = origin._y;
  attr._zorigin = origin._z;
  for (i = 0; i < 3; ++i){
    attr._xaxis[i] = xaxis[i];
    attr._yaxis[i] = yaxis[i];
    attr._zaxis[i] = zaxis[i];
  }
  irtkRealImage image(attr);
//  irtkRealImage image(tgtDim, srcDim, 1, dTgt, dSrc, 1, origin, xaxis, yaxis, zaxis);

  for (j = 0; j < srcDim; ++j){
    for (i = 0; i < tgtDim; ++i){

      x = i; y = j; z = k;
      image.ImageToWorld(x, y, z);

      val = 1000000 * parzenHist.Density(x, y);
      image.Put(i, j, 0, val);
    }
  }

  image.Write("parzen.nii.gz");

  exit(0);

  /////////////////////////////////////////////////
  // Testing not currently used.

//   double density, lhood = 0.0;
//   int trainingCount = 0;
//   int testingCount = 0;

//   for (k = 0; k < zdim; ++k){
//     for (j = 0; j < ydim; ++j){
//       for (i = 0; i < xdim; ++i){
//         x = i;
//         y = j;
//         z = k;
//         targetVal = target->Get(i, j, k);

//         if (targetVal <= padding)
//           continue;

//         target->ImageToWorld(x, y, z);
//         transformation->Transform(x, y, z);
//         source->WorldToImage(x, y, z);

//         if ((x > x1) && (x < x2) &&
//             (y > y1) && (y < y2) &&
//             (z > z1) && (z < z2)){

//           sourceVal = interpolator->EvaluateInside(x, y, z);

//           parzenHist.DelSample(targetVal, sourceVal);
//           density = parzenHist.Density(targetVal, sourceVal);

//           if (density > 0){
//             lhood += log(density);
//           }

//           parzenHist.AddSample(targetVal, sourceVal);
//         }
//       }
//     }
//   }


//   cout << "Training Count: " << trainingCount << endl;
//   cout << "Testing Count:  " << testingCount << endl;
//   cout << "Hist samples:   " << parzenHist.NumberOfSamples() << endl;
//   cout << "Lhood:        " << -1.0 * lhood << endl;


//   float *tgtTesting, *srcTesting;

//   tgtTesting  = new float[xdim * ydim * zdim];
//   srcTesting  = new float[xdim * ydim * zdim];

//   int zerosCount = 0;
//   for (i = 0; i < testingCount; ++i){
//     density = parzenHist.Density(tgtTesting[i], srcTesting[i]);
//     if (density > 0){
//       entropy += density * log(density);
//     } else {
//       ++zerosCount;
//     }
//   }
//   cout << "Zeros count: " << zerosCount << endl;


}


