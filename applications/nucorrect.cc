
#include <sys/types.h>
#include <sys/time.h>

#include <irtkRegistration.h>

// Default filenames
char *gm_name;
char *wm_name;
char *target_name;
char *dofin_name;
char *output_name;

int paddingValue, numberOfLevels;

#define SAMPLES 1000

void fitField(irtkGreyImage &mask, irtkRealImage &logTarget, irtkBSplineFreeFormTransformation *affd);
void fitFieldRandomised(irtkGreyImage &mask, irtkRealImage &logTarget, irtkBSplineFreeFormTransformation *affd);

void usage()
{
  cerr << "nucorrect [image] [gmMask] [wmMask] [output]  <options>\n" << endl;
  cerr << "-ds spacing for lattice used in bias field." << endl;
  cerr << "" << endl;
  cerr << "All input images must have same dimensions." << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{

  // Check command line
  if (argc < 5){
    usage();
  }

  bool ok;
  int  i, j, k, xdim, ydim, zdim;
  double x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis[3], yaxis[3], zaxis[3];
  int voxels, noOfCPs;
  int gmCount, wmCount;
  double x, y, z;
  double temp, curr, avDiff;
  int iterations = 5;

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  gm_name = argv[1];
  argc--;
  argv++;
  wm_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
 
  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkRealImage target(target_name);
  irtkGreyImage output(target_name);
  irtkRealImage logTarget(target_name);
  cout << "done" << endl;

  irtkGreyImage gmMask(gm_name);
  irtkGreyImage wmMask(wm_name);

  // Fix spacing
  dx = 40;
  dy = 40;
  dz = 40;

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)){
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)){
      argc--;
      argv++;
      ok = true;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // Create ffd
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = target.GetX()-1;
  y2 = target.GetY()-1;
  z2 = target.GetZ()-1;

  target.ImageToWorld(x1, y1, z1);
  target.ImageToWorld(x2, y2, z2);
  target.GetOrientation(xaxis, yaxis, zaxis);
    
  irtkBSplineFreeFormTransformation *affd1 = new 
    irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);
  irtkBSplineFreeFormTransformation *affd2 = new 
    irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);
  irtkBSplineFreeFormTransformation *affdAv = new 
    irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);

  xdim = target.GetX();
  ydim = target.GetY();
  zdim = target.GetZ();

  irtkGreyPixel *ptr2gm, *ptr2wm;
  irtkRealPixel *ptr2target, *ptr2logTarget;

  voxels = target.GetNumberOfVoxels();

  // First make sure that every positive mask value refers to a positive
  // intensity in the anatomy, this is because logs will be taken later.

  gmCount = 0;
  ptr2gm = gmMask.GetPointerToVoxels();
  ptr2target = target.GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    if (*ptr2gm > 0){
      if (*ptr2target > 0){
        ++gmCount;
      } else {
        *ptr2gm = 0;
      }
    }
    ++ptr2gm;
    ++ptr2target;
  }

  wmCount = 0;
  ptr2wm = wmMask.GetPointerToVoxels();
  ptr2target = target.GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    if (*ptr2wm > 0){
      if (*ptr2target > 0){
        ++wmCount;
      } else {
        *ptr2wm = 0;
      }
    }
    ++ptr2wm;
    ++ptr2target;
  }

  // Check we have enough data.
  if (gmCount < SAMPLES || wmCount < SAMPLES){
    cerr << "Fewer than " << SAMPLES << " samples to use" << endl;
    exit(1);
  }

  // Take the log of the image where tissue is present.
  ptr2gm = gmMask.GetPointerToVoxels();
  ptr2wm = wmMask.GetPointerToVoxels();
  ptr2target    = target.GetPointerToVoxels();
  ptr2logTarget = logTarget.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*ptr2gm > 0 || *ptr2wm > 0){
      *ptr2logTarget = log(*ptr2target);
    }
    ++ptr2gm;
    ++ptr2wm;
    ++ptr2target;
    ++ptr2logTarget;
  }

//   fitField(gmMask, logTarget, affd1);
//   fitField(wmMask, logTarget, affd2);
  for (int n = 0; n < affdAv->NumberOfDOFs(); ++n){
    affdAv->Put(n, 0.0);
  }

  noOfCPs = affd1->NumberOfDOFs() / 3;
  for (i = 0; i < iterations; ++i){
    fitFieldRandomised(gmMask, logTarget, affd1);
    fitFieldRandomised(wmMask, logTarget, affd2);

    avDiff = 0.0;
    for (int n = 0; n < noOfCPs; ++n){
      temp = affd1->Get(n);
      temp += affd2->Get(n);
      temp *= 0.5;
      temp /= ((double) iterations);

      curr = affdAv->Get(n);
      affdAv->Put(n, temp + curr);
      if (i > 0){
        avDiff += iterations *( fabs((curr / i) - (temp + curr) / (i + 1) ) );
      }
    }
    cout << "Iteration: " << i + 1 << "  Av. Diff of CP values for field " << avDiff / (double) noOfCPs << endl; 
  }

  // The values of the displacement field are an estimate for the log of
  // the bias field.  I.e. observed = true * field which implies true =
  // observed / field, i.e. true = observed * exp (-1.0 * logfield).

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){

        x = i;
        y = j;
        z = k;
        target.ImageToWorld(x, y, z);
        affdAv->LocalDisplacement(x, y, z);

        temp = target.Get(i, j, k);
        temp = temp * exp(-1.0 * x);
        output.Put(i, j, k, round(temp));

      }
    }
  }

  output.Write(output_name);

}

/////////////////////////////////////////////////////////////////////////////////////////////////

void fitFieldRandomised(irtkGreyImage &mask, irtkRealImage &logTarget, irtkBSplineFreeFormTransformation *affd){

//   cout << endl;
//   cout << "In fitFieldRandomised function." << endl;

  irtkGreyPixel *ptr2mask;
  irtkRealPixel *ptr2logTarget;
  int maskCount, n, voxels, index;
  int i, j, k, xdim, ydim, zdim, noOfDofs;
  double x, y, z, meanVal, fraction;

  double xlocs[SAMPLES], ylocs[SAMPLES], zlocs[SAMPLES];
  double vals[SAMPLES], dummya[SAMPLES], dummyb[SAMPLES];
  // dummy variables for a FFD approximate call.
  for (i = 0; i < SAMPLES; ++i){
    dummya[i] = dummyb[i] = 0.0;
  }

  timeval tv;
  gettimeofday(&tv, NULL); 
  long init = tv.tv_usec;
  long temp = -1*init;
  (void) ran2(&temp);

  // blah = ran2(&init); // gets a random number.

  xdim = logTarget.GetX();
  ydim = logTarget.GetY();
  zdim = logTarget.GetZ();

  voxels = mask.GetNumberOfVoxels();

  // How many voxels in the mask.
  maskCount = 0;
  ptr2mask = mask.GetPointerToVoxels();
  for (n = 0; n < voxels; ++n){
    if (*ptr2mask > 0){
      ++maskCount;
    }
    ++ptr2mask;
  }

  // What fraction of the mask ends up selected?
  fraction = ((double) SAMPLES) / maskCount;

  // Some reporting.
//   cout << "maskCount " << maskCount << endl;
//   cout << "fraction : " << fraction << endl;

  // Subsample the data and store it for the FFD approximation.

  index     = 0;
  while (index < SAMPLES){

    ptr2mask = mask.GetPointerToVoxels();
    ptr2logTarget = logTarget.GetPointerToVoxels();

    for (n = 0; n < voxels; ++n){

      if (*ptr2mask > 0){

        // randomly decide whether or not to include sample
        if (ran2(&init) <= fraction){
          // include this in the approximation.
          vals[index] = *ptr2logTarget;
          // Get the world coordinates
          i = n % xdim;
          j = (n % (xdim * ydim)) / xdim ;
          k = n / (xdim * ydim);
          x = i;
          y = j;
          z = k;

          logTarget.ImageToWorld(x, y, z);
          xlocs[index] = x;
          ylocs[index] = y;
          zlocs[index] = z;

          ++index;

          if (index >= SAMPLES){
            break;
          }
        }
      }

      ++ptr2mask;
      ++ptr2logTarget;
    }
  }


  // Mean centre the data before passing to the FFD approximate function.
  // This is so that the FFD is not required to approximate large values
  // int ehregion of the mask.
  meanVal = 0;
  for (n = 0; n < SAMPLES; ++n){
    meanVal += vals[n];
  }
  meanVal /= SAMPLES;
//   cout << "Mean val pre approx: " << meanVal  << "  Centering ..." << endl;
  for (n = 0; n < SAMPLES; ++n){
    vals[n] -= meanVal;
  }

  affd->Approximate(xlocs, ylocs, zlocs, vals, dummya, dummyb, SAMPLES);

  // Mean centre the ffd's x displacements after approximation.
  noOfDofs = affd->NumberOfDOFs();
  meanVal = 0;

  for (i = 0; i < noOfDofs / 3; ++i){
    meanVal += affd->Get(i);
  }

  meanVal = meanVal / (noOfDofs / 3);
//   cout << "Mean FFD x disp after approx : " << meanVal << endl;

  for (i = 0; i < noOfDofs / 3; ++i){
    affd->Put(i, affd->Get(i) - meanVal);
  }


  double error = 0;

  for (n = 0; n < SAMPLES; ++n){
    x = xlocs[n];
    y = ylocs[n];
    z = zlocs[n];
    affd->LocalDisplacement(x, y, z);
    error += (x - vals[n]) * (x - vals[n]);
  }

  error /= SAMPLES;

//   cout << "MSE : " << error << endl;
//   cout << "RMS : " << sqrt(error) << endl;
//   cout << endl;

}

void fitField(irtkGreyImage &mask, irtkRealImage &logTarget, irtkBSplineFreeFormTransformation *affd){

  cout << endl;
  cout << "In fitField function." << endl;

  irtkGreyPixel *ptr2mask;
  irtkRealPixel *ptr2logTarget;
  int maskCount, n, voxels, step, index;
  int i, j, k, xdim, ydim, zdim, noOfDofs;
  double x, y, z, meanVal;

  double xlocs[SAMPLES], ylocs[SAMPLES], zlocs[SAMPLES];
  double vals[SAMPLES], dummya[SAMPLES], dummyb[SAMPLES];
  // dummy variables for a FFD approximate call.
  for (i = 0; i < SAMPLES; ++i){
    dummya[i] = dummyb[i] = 0.0;
  }

  xdim = logTarget.GetX();
  ydim = logTarget.GetY();
  zdim = logTarget.GetZ();

  voxels = mask.GetNumberOfVoxels();

  // How many voxels in the mask.
  maskCount = 0;
  ptr2mask = mask.GetPointerToVoxels();
  for (n = 0; n < voxels; ++n){
    if (*ptr2mask > 0){
      ++maskCount;
    }
    ++ptr2mask;
  }

  // Step to determine how to downsample the data.
  step = maskCount / SAMPLES;

  // Subsample the data and store it for the FFD approximation.
  ptr2mask = mask.GetPointerToVoxels();
  ptr2logTarget = logTarget.GetPointerToVoxels();

  maskCount = 0;
  index     = 0;

  for (n = 0; n < voxels; ++n){

    // Ignore the first voxel.
    if (*ptr2mask > 0 && maskCount == 0){
      ++maskCount;
      continue;
    }

    if (*ptr2mask > 0 && maskCount % step == 0){
      ++maskCount;

      // include this in the approximation.
      vals[index] = *ptr2logTarget;

      // Get the world coordinates
      i = n % xdim;
      j = (n % (xdim * ydim)) / xdim ;
      k = n / (xdim * ydim);

      x = i;
      y = j;
      z = k;

      logTarget.ImageToWorld(x, y, z);

      xlocs[index] = x;
      ylocs[index] = y;
      zlocs[index] = z;

      ++index;

    } else if (*ptr2mask > 0){
      // an ignored voxel
      ++maskCount;
    }

    ++ptr2mask;
    ++ptr2logTarget;
  }

  // Some reporting.
  cout << "maskCount " << maskCount << endl;
  cout << "sub-sampled data in steps of " << step << endl;

  // Mean centre the data before passing to the FFD approximate function.
  // This is so that the FFD is not required to approximate large values
  // int ehregion of the mask.
  meanVal = 0;
  for (n = 0; n < SAMPLES; ++n){
    meanVal += vals[n];
  }
  meanVal /= SAMPLES;
  cout << "Mean val pre approx: " << meanVal  << "  Centering ..." << endl;
  for (n = 0; n < SAMPLES; ++n){
    vals[n] -= meanVal;
  }

  affd->Approximate(xlocs, ylocs, zlocs, vals, dummya, dummyb, SAMPLES);

  // Mean centre the ffd's x displacements after approximation.
  noOfDofs = affd->NumberOfDOFs();
  meanVal = 0;

  for (i = 0; i < noOfDofs / 3; ++i){
    meanVal += affd->Get(i);
  }

  meanVal = meanVal / (noOfDofs / 3);
  cout << "Mean FFD x disp after approx : " << meanVal << endl;

  for (i = 0; i < noOfDofs / 3; ++i){
    affd->Put(i, affd->Get(i) - meanVal);
  }


  double error = 0;

  for (n = 0; n < SAMPLES; ++n){
    x = xlocs[n];
    y = ylocs[n];
    z = zlocs[n];
    affd->LocalDisplacement(x, y, z);
    error += (x - vals[n]) * (x - vals[n]);
  }

  error /= SAMPLES;

  cout << "MSE : " << error << endl;
  cout << "RMS : " << sqrt(error) << endl;
  cout << endl;

}




//   if (0){
//   // Step to determine how to downsample the data.
//   gmStep = gmCount / SAMPLES;
//   wmStep = wmCount / SAMPLES;

//   // Subsample the data and store it for the FFD approximation.
//   ptr2gm = gmMask.GetPointerToVoxels();
//   ptr2logTarget = logTarget.GetPointerToVoxels();
//   gmCount = 0;

//   for (n = 0; n < voxels; ++n){

//     // Ignore the first voxel.
//     if (*ptr2gm > 0 && gmCount == 0){
//       ++gmCount;
//       continue;
//     }

//     if (*ptr2gm > 0 && gmCount % gmStep == 0){
//       ++gmCount;

//       // include this in the approximation.
//       vals[index] = *ptr2logTarget;

//       // Get the world coordinates
//       i = n % xdim;
//       j = (n % (xdim * ydim)) / xdim ;
//       k = n / (xdim * ydim);

//       x = i;
//       y = j;
//       z = k;

//       logTarget.ImageToWorld(x, y, z);

//       xlocs[index] = x;
//       ylocs[index] = y;
//       zlocs[index] = z;

//       ++index;

//     } else if (*ptr2gm > 0){
//       // an ignored voxel
//       ++gmCount;
//     }

//     ++ptr2gm;
//     ++ptr2logTarget;
//   }

//   // Some reporting.
//   cout << "gmCount " << gmCount << " index " << index << endl;

//   // Mean centre the data before passing to the FFD approximate function.
//   // This is so that the FFD is not required to approximate large values
//   // int ehregion of the mask.
//   meanVal = 0;
//   for (n = 0; n < SAMPLES; ++n){
//     meanVal += vals[n];
//   }
//   meanVal /= SAMPLES;
//   cout << "Mean val pre approx: " << meanVal  << "  Centering ..." << endl;
//   for (n = 0; n < SAMPLES; ++n){
//     vals[n] -= meanVal;
//   }

//   affd->Approximate(xlocs, ylocs, zlocs, vals, dummya, dummyb, SAMPLES);

//   // Mean centre the ffd's x displacements after approximation.
//   noOfDofs = affd->NumberOfDOFs();
//   meanVal = 0;

//   for (i = 0; i < noOfDofs / 3; ++i){
//     meanVal += affd->Get(i);
//   }

//   meanVal = meanVal / (noOfDofs / 3);
//   cout << "Mean FFD x disp after approx : " << meanVal << endl;

//   for (i = 0; i < noOfDofs / 3; ++i){
//     affd->Put(i, affd->Get(i) - meanVal);
//   }

//   } // if (0)



//   exit(0);

//   meanVal = 0;
//   double ssq = 0;
//   for (i = 0; i < affd->GetX(); ++i){
//     for (j = 0; j < affd->GetY(); ++j){
//       for (k = 0; k < affd->GetZ(); ++k){
//         affd->Get(i, j, k, x, y, z);
//         meanVal += x;
//         ssq += x * x;
//       }
//     }
//   }
//   n = affd->GetX() * affd->GetY() * affd->GetZ();

//   meanVal /= n;
//   cout << "meanVal " << meanVal << endl;
//   cout << "var     " << (ssq / ((double) n)) - (meanVal * meanVal) << endl;
//   mffd->irtkTransformation::Write("temp.dof");



//   irtkBiasCorrection biasCorrection;
  
//   // Set input and output for the biasCorrection filter
//   biasCorrection.SetInput(&target, &target);
//   biasCorrection.SetOutput(mffd);

//   biasCorrection.SetGreyMask(&gmMask);
//   biasCorrection.SetWhiteMask(&wmMask);

//   if (parameter_name != NULL){
//     biasCorrection.Read(parameter_name);
//   }

//   // Run biasCorrection filter
//   biasCorrection.Run();

//   // write a corrected version of the target.
//   biasCorrection.CorrectTarget();

//   target.Write(output_name);
