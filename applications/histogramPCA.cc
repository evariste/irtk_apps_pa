#include <irtkImage.h>
#include <irtkHistogram.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>

char *input_name1 = NULL, *output_name = NULL;
char *input_name2 = NULL;
int _stepSize = 1;

void usage(){

  cout << "histogramPCA [imageA] [imageB] [radius]" << endl;
  cout << "-TpA: padding A. -TpB: padding B." << endl;
  cout << "-output: name for output image containing histogram, if required." << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;
  exit(0);
}

// Displace range of unpadded voxels to start at zero.
// Padded voxels set to zero.
void preProcess(irtkGreyImage *image, irtkGreyPixel pad){
  int i, voxels;
  irtkGreyPixel *pPix;
  irtkGreyPixel min;

  voxels = image->GetNumberOfVoxels();
  pPix   = image->GetPointerToVoxels();
  min    = MAX_GREY;

  for (i = 0; i < voxels; ++i){
    if (*pPix > pad){
      if (min > *pPix)
        min = *pPix;
    }
    ++pPix;
  }

  pPix = image->GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    if (*pPix > pad){
      *pPix -= min;
    } else {
      *pPix = 0;
    }
    ++pPix;
  }
}

// Also generates the mask.
int getWindowCount(irtkGreyImage *image1, irtkGreyImage *image2, 
                   int radius, irtkGreyImage *mask){
  int i, j, k;
  int a, b, c;
  int count = 0;
  int minVal1, maxVal1, minVal2, maxVal2;

  irtkGreyPixel *pCentre1, *pCentre2, *pTemp1, *pTemp2;
  irtkGreyPixel *pMask;
  irtkGreyPixel high1, low1;
  irtkGreyPixel high2, low2;

  int xmax = image1->GetX() - 1 - radius;
  int ymax = image1->GetY() - 1 - radius;
  int zmax = image1->GetZ() - 1 - radius;
  int xoff = 1;
  int yoff = image1->GetX();
  int zoff = image1->GetX() * image1->GetY();

  image1->GetMinMax(&low1, &high1);
  image2->GetMinMax(&low2, &high2);

  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre1 = image1->GetPointerToVoxels(radius, j, k);
      pCentre2 = image2->GetPointerToVoxels(radius, j, k);
      pMask    = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){

        minVal1 = high1;
        maxVal1 = low1;
        minVal2 = high2;
        maxVal2 = low2;

        for (c = -1 * radius; c <= radius; ++c){
          for (b = -1 * radius; b <= radius; ++b){
            for (a = -1 * radius; a <= radius; ++a){
              pTemp1 = pCentre1 + a * xoff + b * yoff + c * zoff;
              pTemp2 = pCentre2 + a * xoff + b * yoff + c * zoff;

              if (minVal1 > *pTemp1)
                minVal1 = *pTemp1;
              if (maxVal1 < *pTemp1)
                maxVal1 = *pTemp1;
              if (minVal2 > *pTemp2)
                minVal2 = *pTemp2;
              if (maxVal2 < *pTemp2)
                maxVal2 = *pTemp2;
            }
          }
        }

        if ((minVal1 > 0 && maxVal1 > minVal1) ||
            (minVal2 > 0 && maxVal2 > minVal2) ){
          ++count;
          *pMask = 1;
        } else {
          *pMask = 0;
        }

        pCentre1 += _stepSize;
        pCentre2 += _stepSize;
        pMask += _stepSize;
      }
    }
  }

  return count;
}
//////////////////////////////////////////////////////////////////
// Gets the mean window across a single image
void getMeanVector(irtkGreyImage *image,
                   irtkGreyImage *mask,
                   int radius, int dims,
                   irtkVector &meanWindow){
  int i, j, k;
  int a, b, c;
  int index, count = 0;

  irtkGreyPixel *pCentre, *pTemp, *pMask;

  int xmax = image->GetX() - 1 - radius;
  int ymax = image->GetY() - 1 - radius;
  int zmax = image->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = image->GetX();
  int zoff = image->GetX() * image->GetY();

  // Initialise mean.
  for (i = 0; i < dims; ++i)
    meanWindow(i) = 0.0;
  // Main loop.  
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre = image->GetPointerToVoxels(radius, j, k);
      pMask    = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
        if(*pMask > 0){
          index = 0;
          //Loop over window in first image.
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp = pCentre + a * xoff + b * yoff + c * zoff;
                meanWindow(index) += *pTemp;
                ++index;
              }
            }
          }
          count++;
        }
        pCentre += _stepSize;
        pMask += _stepSize;
      }
    }
  }
  meanWindow = meanWindow * (1.0 / count);
}
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// This routine avoids the overhead of storing the data for all windows in
// a single (huge) data matrix.
void getCovMatrix(irtkMatrix &covMatrix,
                  irtkGreyImage *image,
                  irtkGreyImage *mask,
                  int radius, int dims,
                  irtkVector &meanWindow){
  int i, j, k;
  int a, b, c, d1, d2;
  int count = 0;
  int index;
  irtkGreyPixel *pCentre, *pTemp, *pMask;
  irtkMatrix temp;
  irtkVector currVector;
  temp.Initialize(dims, dims);
  currVector.Initialize(dims);
  // Limits.
  int xmax = image->GetX() - 1 - radius;
  int ymax = image->GetY() - 1 - radius;
  int zmax = image->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = image->GetX();
  int zoff = image->GetX() * image->GetY();

  // Main loop.
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre = image->GetPointerToVoxels(radius, j, k);
      pMask   = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
        if(*pMask > 0){
          //Valid local windows.
          currVector *= 0.0;
          index = 0;
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp = pCentre + a * xoff + b * yoff + c * zoff;
                currVector(index) = *pTemp;
                ++index;
              }
            }
          }
          ++count;
          // Find x - mu
          currVector -= meanWindow;

          // add (x - mu)^T * (x - mu) to the accumulating cov matrix.
          for (d1 = 0; d1 < dims; ++d1){
            for (d2 = d1; d2 < dims; ++d2){
              temp(d1, d2) = temp(d2, d1) = currVector(d1) * currVector(d2);
            }
          }
          covMatrix += temp;
        } // End of valid local windows section.

        pCentre += _stepSize;
        pMask += _stepSize;
      }
    }
  }
  // Normalise the cov matrix.
  covMatrix = covMatrix * (1.0 / (count - 1));
}
/////////////////////////////////////////////////////////////////


bool equalImageSizes(irtkGreyImage *first, irtkGreyImage *second){
  return ( first->GetX() == second->GetX() ) &&
         ( first->GetY() == second->GetY() ) &&
         ( first->GetZ() == second->GetZ() );
}

int countMainModes(irtkEigenAnalysis &ea, int dims, float minPercentage){

  int i, mainModeCount = 0;
  float fTotalVar   = 0.0;
  //  float fCumulated = 0.0;

  // Find sum of eigenvalues.
  for (i = 0; i < dims; i++){
    fTotalVar += ea.Eigenvalue(i);
  }

  // Identify the main modes.
  for (i = 0; i < dims; i++){

    if (100 * ea.Eigenvalue(i) / fTotalVar < minPercentage)
      break;

    mainModeCount++;
//     cout<< "Mode:" << i+1 << " Eigenvalue:" << ea.Eigenvalue(i) << endl;
//     fCumulated += 100 * ea.Eigenvalue(i) / fTotalVar;
//     cout << " Mode [" << i+1 << "] explains " 
//          << 100 * ea.Eigenvalue(i) / fTotalVar
//          << " % (" 
//          << fCumulated 
//          << " %) of shape variance." << endl;
  }

  cout  << "mainModeCount :" << mainModeCount << endl;
  return mainModeCount;
}

void getProjectedData (irtkMatrix &projData, irtkMatrix &mainModeVects,
                       irtkGreyImage *image, int mainModeCount,
                       int windowCount, int dims, int radius, irtkGreyImage *mask,
                       irtkVector &meanWindow
){
  int i, j, k;
  int a, b, c, mode, dim;
  int count = 0;
  int index;
  irtkGreyPixel *pCentre, *pTemp, *pMask;
  irtkVector currVector;
  currVector.Initialize(dims);

  // Limits.
  int xmax = image->GetX() - 1 - radius;
  int ymax = image->GetY() - 1 - radius;
  int zmax = image->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = image->GetX();
  int zoff = image->GetX() * image->GetY();

  // Main loop.
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre = image->GetPointerToVoxels(radius, j, k);
      pMask   = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
        if(*pMask > 0){
          //Valid local windows.
          currVector *= 0.0;
          index = 0;
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp = pCentre + a * xoff + b * yoff + c * zoff;
                currVector(index) = *pTemp;
                ++index;
              }
            }
          }

          // Find x - mu
          currVector -= meanWindow;

          for (mode = 0; mode <  mainModeCount; ++mode){
            projData(mode, count) = 0.0;
            for (dim = 0; dim < dims; ++dim){
              projData(mode, count) += mainModeVects(dim, mode) * currVector(dim);
            }
          }

          ++count;
        } // End of valid local windows section.

        pCentre += _stepSize;
        pMask += _stepSize;
      }
    }
  }
}

void normaliseProjData (irtkMatrix &projData, irtkVector &eigenvals, 
                        int mainModeCount, int dataItems){
  // projData contains the components of the original data w.r.t. the main
  // eigen modes.
  // Divide each component by the corresponding eigenvalue (variance).
  // This is done to get the Mahalanobis distance later.
  int i, j;
  float variance;

  for (i = 0; i < mainModeCount; ++i){
    variance = eigenvals(i);
    for (j = 0; j < dataItems;     ++j){
      projData(i, j) = projData(i, j) / variance;
    }
  }
}

void getDistances (irtkVector &dists, irtkMatrix &projData, int mainModeCount, int dataItems){
  int i, j;

  for (j = 0; j < dataItems;     ++j){
    dists(j) = 0.0;
    for (i = 0; i < mainModeCount; ++i){
      dists(j) += projData(i, j) * projData(i, j);
    }
    dists(j) = sqrt(dists(j));
  }
}

float getMaxDist(irtkVector &dists, int dataItems){
  int i;
  float max = FLT_MIN;
  for (i = 0; i < dataItems; ++i){
    if (max < dists(i)){
      max = dists(i);
    }
  }
  return max;
}

int main(int argc, char **argv)
{
  irtkMatrix covMatrix1, covMatrix2;
  irtkMatrix mainModeVects1, mainModeVects2;
  irtkMatrix projData1, projData2;
  irtkVector meanWindow1, meanWindow2;
  irtkVector eigenvals1, eigenvals2;
  irtkVector dists1, dists2;

  int radius = 1, i, j;
  bool ok;
  irtkGreyImage *input1, *input2, *mask;
  irtkGreyPixel padding1 = MIN_GREY, padding2 = MIN_GREY;
  irtkRealImage *output;

  int   windowCount, dims;
  int   mainModeCount1 = 0;
  int   mainModeCount2 = 0;
  float minPercentage = 1.0;
  int   bins = 100;
  float maxDist1, maxDist2;
  bool overlapping = false;

  // Parse arguments.
  if (argc < 4){
    usage();
  }

  input_name1 = argv[1];
  argv++; argc--;
  input_name2 = argv[1];
  argv++; argc--;
  radius      = atoi(argv[1]);
  argv++; argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-TpA") == 0)){
      argc--;      argv++;
      padding1 = atoi(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-TpB") == 0)){
      argc--;      argv++;
      padding2 = atoi(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
      output_name = argv[1];
      argv++; argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-overlap") == 0)){
      overlapping = true;
      argv++; argc--;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Decide if the windows are overlapping or tiled (default).
  if (overlapping == true){
    _stepSize = 1;
  } else {
    _stepSize = radius;
  }

  // Read and pre-process images.
  input1 = new irtkGreyImage(input_name1);
  input2 = new irtkGreyImage(input_name2);
  mask   = new irtkGreyImage(input_name1);

  if (!equalImageSizes(input1, input2)){
    cout << "Input images: sizes don't match" << endl;
    exit(1);
  }

  preProcess(input1, padding1);
  preProcess(input2, padding2);

  // Count the number of valid windows and generate the mask.
  windowCount = getWindowCount(input1, input2, radius, mask);

  // dims is the number of voxels in a window.
  dims = (int) pow(2 * radius + 1, 3.0);

  meanWindow1.Initialize(dims);
  meanWindow1 *= 0.0;
  getMeanVector(input1, mask, radius, dims, meanWindow1);
  covMatrix1.Initialize(dims, dims);
  getCovMatrix(covMatrix1, input1, mask, radius, dims, meanWindow1);
  // Compute the eigen vals and vecs for the covariance matrix.
  irtkEigenAnalysis ea1(dims);
  for (i = 0; i < dims; ++i){
    for(j = 0; j < dims; ++j){
      ea1.Matrix(i, j) = covMatrix1(i, j);  
    }
  }
  ea1.DecrSortEigenStuff();

  // How many modes contribute at lease the min percentage of the variance?
  mainModeCount1 = countMainModes(ea1, dims, minPercentage);
  // Copy the main modes
  mainModeVects1.Initialize(dims, mainModeCount1);
  mainModeVects1 *= 0.0;
  for (i = 0; i < dims ;          ++i){
    for (j = 0; j < mainModeCount1; ++j){
      mainModeVects1(i, j) = ea1.Eigenvector(i, j);
    }
  }
  // Copy the main eigenvalues.
  eigenvals1.Initialize(mainModeCount1);
  for (i = 0; i < mainModeCount1 ; i++){
    eigenvals1(i) = ea1.Eigenvalue(i);
  }

  // The data is projected onto a basis consisting of the main mode
  // vectors. The size is mainModeCount x the number of data items.
  projData1.Initialize(mainModeCount1, windowCount);
  getProjectedData(projData1, mainModeVects1, input1,
                   mainModeCount1, windowCount, dims, radius, mask, meanWindow1);

  normaliseProjData(projData1, eigenvals1, mainModeCount1, windowCount);

  dists1.Initialize(windowCount);
  getDistances(dists1, projData1, mainModeCount1, windowCount);
  // Free memory
  projData1.Initialize(1, 1);

  maxDist1 = getMaxDist(dists1, windowCount);

  //////////////////////////////////////////////////////////////////////
  /// image 2.
  /////////////////////////////////////////////////////////////////////

  // Initialise data structures.
  meanWindow2.Initialize(dims);
  getMeanVector(input2, mask, radius, dims, meanWindow2);
  covMatrix2.Initialize(dims, dims);
  getCovMatrix(covMatrix2, input2, mask, radius, dims, meanWindow2);
  // Compute the eigen vals and vecs for the covariance matrices.
  irtkEigenAnalysis ea2(dims);
  for (i = 0; i < dims; ++i){
    for(j = 0; j < dims; ++j){
      ea2.Matrix(i, j) = covMatrix2(i, j);
    }
  }
  ea2.DecrSortEigenStuff();

  mainModeCount2 = countMainModes(ea2, dims, minPercentage);
  // Copy the main modes
  mainModeVects2.Initialize(dims, mainModeCount2);
  mainModeVects2 *= 0.0;
  for (i = 0; i < dims ;          ++i){
    for (j = 0; j < mainModeCount2; ++j){
      mainModeVects2(i, j) = ea2.Eigenvector(i, j);
    }
  }

  // Copy the main eigenvalues.
  eigenvals2.Initialize(mainModeCount2);
  for (i = 0; i < mainModeCount2 ; i++){
    eigenvals2(i) = ea2.Eigenvalue(i);
  }

  // The data is projected onto a basis consisting of the main mode
  // vectors. The size is mainModeCount x the number of data items.
  projData2.Initialize(mainModeCount2, windowCount);
  getProjectedData(projData2, mainModeVects2, input2,
                   mainModeCount2, windowCount, dims, radius, mask, meanWindow2);

  normaliseProjData(projData2, eigenvals2, mainModeCount2, windowCount);

  dists2.Initialize(windowCount);
  getDistances(dists2, projData2, mainModeCount2, windowCount);
  // Free memory.
  projData2.Initialize(1, 1);

  maxDist2 = getMaxDist(dists2, windowCount);

  // Now the main bit.
  irtkHistogram_2D histogram(bins, bins);

  histogram.PutMin(0,  0);
  histogram.PutMax(maxDist1, maxDist2);
  for (i = 0; i < windowCount; ++i){
    histogram.AddSample(dists1(i), dists2(i));
  }

  output = new irtkRealImage(bins, bins, 1, 1.0, 1.0, 1.0);
  for (j = 0; j < bins; ++j){
    for (i = 0; i < bins; ++i){
      output->Put(i, j, 0, histogram(i, j));
    }
  }

  if (output_name != NULL){
    output->Write(output_name);
  }

  cout << "MI: " << histogram.MutualInformation() << endl;
  cout << "NMI: " << histogram.NormalizedMutualInformation() << endl;

}


// // Populates data from two images simultaneously.
// // Data matrix and mean vector assumed to have been initialised
// // to correct dimensions.
// int getDataAndMean (irtkGreyImage *image1, irtkGreyImage *image2,
//                     irtkGreyImage *mask,
//                     int radius, int dims,
//                     irtkMatrix &dataMatrix, irtkVector &meanWindow){
//   int i, j, k;
//   int a, b, c;
//   int count = 0;
//   int index;

//   irtkGreyPixel *pCentre1, *pTemp1, *pCentre2, *pTemp2, *pMask;

//   int xmax = image1->GetX() - 1 - radius;
//   int ymax = image1->GetY() - 1 - radius;
//   int zmax = image1->GetZ() - 1 - radius;
//   // Offsets.
//   int xoff = 1;
//   int yoff = image1->GetX();
//   int zoff = image1->GetX() * image1->GetY();

//   // Initialise mean.
//   for (i = 0; i < dims; ++i){
//     meanWindow(i) = 0.0;
//   }

//   for (k = radius; k < zmax; ++k){
//     for (j = radius; j < ymax; ++j){

//       pCentre1 = image1->GetPointerToVoxels(radius, j, k);
//       pCentre2 = image2->GetPointerToVoxels(radius, j, k);
//       pMask    = mask->GetPointerToVoxels(radius, j, k);

//       for (i = radius; i < xmax; ++i){

//         if(*pMask > 0){

//           index = 0;
//           for (c = -1 * radius; c <= radius; ++c){
//             for (b = -1 * radius; b <= radius; ++b){
//               for (a = -1 * radius; a <= radius; ++a){
//                 pTemp1 = pCentre1 + a * xoff + b * yoff + c * zoff;

//                 dataMatrix(index, count) = *pTemp1;
//                 meanWindow(index) += *pTemp1;
//                 ++index;
//               }
//             }
//           }

//           for (c = -1 * radius; c <= radius; ++c){
//             for (b = -1 * radius; b <= radius; ++b){
//               for (a = -1 * radius; a <= radius; ++a){
//                 pTemp2 = pCentre2 + a * xoff + b * yoff + c * zoff;

//                 dataMatrix(index, count) = *pTemp2;
//                 meanWindow(index) += *pTemp2;
//                 ++index;
//               }
//             }
//           }

//           ++count;
//         } // if *pMask > 0

//         ++pCentre1;
//         ++pCentre2;
//         ++pMask;
//       }

//     }
//   }

//   meanWindow = meanWindow * (1.0 / count);

//   // Subtract the mean from the data.
//   for (i = 0; i < count; ++i){
//     for (j = 0; j < dims;  ++j){
//       dataMatrix(j,i) -= meanWindow(j);
//     }
//   }

//   return count;
// }


// // Get data from a single image into the data matrix.
// // Find mean vector.
// int getDataAndMean (irtkGreyImage *image,   irtkGreyImage *mask,
//                     int radius,            int dims,
//                     irtkMatrix &dataMatrix, irtkVector &meanWindow){

//   // Pre: dataMatrix and meanWindow have space allocated.

//   int i, j, k;
//   int a, b, c;
//   int count = 0;
//   int index;

//   irtkGreyPixel *pCentre, *pTemp, *pMask;

//   int xmax = image->GetX() - 1 - radius;
//   int ymax = image->GetY() - 1 - radius;
//   int zmax = image->GetZ() - 1 - radius;

//   int xoff = 1;
//   int yoff = image->GetX();
//   int zoff = image->GetX() * image->GetY();

//   // Init mean.
//   for (i = 0; i < dims; ++i){
//     meanWindow(i) = 0.0;
//   }

//   for (k = radius; k < zmax; ++k){
//     for (j = radius; j < ymax; ++j){

//       pCentre = image->GetPointerToVoxels(radius, j, k);
//       pMask    = mask->GetPointerToVoxels(radius, j, k);

//       for (i = radius; i < xmax; ++i){

//         if(*pMask > 0){

//           index = 0;
//           for (c = -1 * radius; c <= radius; ++c){
//             for (b = -1 * radius; b <= radius; ++b){
//               for (a = -1 * radius; a <= radius; ++a){
//                 pTemp = pCentre + a * xoff + b * yoff + c * zoff;

//                 dataMatrix(index, count) = *pTemp;
//                 meanWindow(index) += *pTemp;
//                 ++index;
//               }
//             }
//           }

//           ++count;
//         }

//         ++pCentre;
//         ++pMask;
//       }

//     }
//   }

//   meanWindow = meanWindow * (1.0 / count);

//   // Subtract the mean from the data.
//   for (i = 0; i < count; ++i){
//     for (j = 0; j < dims;  ++j){
//       dataMatrix(j,i) -= meanWindow(j);
//     }
//   }

//   return count;
// }



// void getCovMatrix(irtkMatrix &dataMatrix, irtkMatrix &T, int dims, int windowCount){

//   // Pre: data matrix has had mean subtracted.
//   int i, j, k;

//   for (i = 0; i < dims; ++i){
//     for (j = 0; j < dims; ++j){
//       T(i, j) = 0.0;
//       for (k = 0; k < windowCount; ++k)
//         T(i, j) += dataMatrix(i, k) * dataMatrix(j, k);
//     }
//   }
//   T = T * (1.0 / (windowCount - 1));

// }
