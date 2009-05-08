#include <irtkImage.h>
#include <irtkEigenAnalysis.h>

char *input_name = NULL, *output_name = NULL;
int _stepSize = 1;

void usage(){

  cout << "pcaDistanceMap [image] [output] [radius]" << endl;
  cout << "-Tp: padding." << endl;
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
  cout << "preproc done" << endl;
}

void getMask(irtkGreyImage *image, int radius, irtkGreyImage *mask, int &windowCount){
  int i, j, k;
  int a, b, c;
  int minVal, maxVal;
  int count = 0, voxels;

  irtkGreyPixel *pCentre, *pTemp;
  irtkGreyPixel *pMask;
  irtkGreyPixel high, low;

  int xmax = image->GetX() - 1 - radius;
  int ymax = image->GetY() - 1 - radius;
  int zmax = image->GetZ() - 1 - radius;
  int xoff = 1;
  int yoff = image->GetX();
  int zoff = image->GetX() * image->GetY();

  image->GetMinMax(&low, &high);
  voxels = image->GetNumberOfVoxels();

  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre = image->GetPointerToVoxels(radius, j, k);
      pMask   = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){

        minVal = high;
        maxVal = low;

        for (c = -1 * radius; c <= radius; ++c){
          for (b = -1 * radius; b <= radius; ++b){
            for (a = -1 * radius; a <= radius; ++a){
              pTemp = pCentre + a * xoff + b * yoff + c * zoff;

              if (minVal > *pTemp)
                minVal = *pTemp;
              if (maxVal < *pTemp)
                maxVal = *pTemp;
            }
          }
        }

        if ( minVal > 0 && maxVal > minVal ){
          *pMask = 1;
          ++count;
        } else {
          *pMask = 0;
        }

        ++pCentre;
        ++pMask;
      }
    }
  }
  cout << "Mask percentage :" << 100 * ((double) count / voxels) << endl;
  windowCount = count;
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
      pMask = mask->GetPointerToVoxels(radius, j, k);

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
        ++pCentre;
        ++pMask;
      }
    }
  }
  meanWindow = meanWindow * (1.0 / count);
  cout << "got mean." << endl;
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
      pMask = mask->GetPointerToVoxels(radius, j, k);

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

        ++pCentre;
        ++pMask;
      }
    }
  }
  // Normalise the cov matrix.
  covMatrix = covMatrix * (1.0 / (count - 1));

  cout << "Got cov matrix." << endl;
}
/////////////////////////////////////////////////////////////////

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
                       irtkGreyImage *image, 
                       irtkGreyImage *mask,
                       int mainModeCount,
                       int dims, int radius,
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
      pMask = mask->GetPointerToVoxels(radius, j, k);

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

        ++pCentre;
        ++pMask;
      }
    }
  }
  cout << "data projected onto main mode basis" << endl;
}
// Given the data projected onto a basis consisting of the main mode
// vectors, use the components relative to this basis to recombine the main
// modes into an image.
void redrawImage (irtkMatrix &projData, irtkMatrix &mainModeVects,
                  irtkRealImage *output,
                  irtkGreyImage *mask,
                  int mainModeCount,
                  int dims, int radius,
                  irtkVector &meanWindow ){
  int i, j, k;
  int a, b, c, mode, dim;
  int count = 0;
  int index;
  irtkRealPixel *pOut, *pTemp;
  irtkGreyPixel *pCounts;
  irtkGreyPixel *pMask, *pTemp2;
  irtkVector currVector;
  currVector.Initialize(dims);

  irtkGreyImage *visitCounts = new irtkGreyImage(input_name);
  pCounts = visitCounts->GetPointerToVoxels();
  int voxels = visitCounts->GetNumberOfVoxels();
  for (i = 0; i < voxels; ++i){
    *pCounts = 0;
    ++pCounts;
  }

  // Limits.
  int xmax = output->GetX() - 1 - radius;
  int ymax = output->GetY() - 1 - radius;
  int zmax = output->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = output->GetX();
  int zoff = output->GetX() * output->GetY();
  // Main loop.
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pOut  = output->GetPointerToVoxels(radius, j, k);
      pMask = mask->GetPointerToVoxels(radius, j, k);
      pCounts = visitCounts->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
        if(*pMask > 0){
          //Valid local windows.
          currVector = meanWindow;
          for (mode = 0; mode <  mainModeCount; ++mode){
            for (dim = 0; dim < dims; ++dim){
              currVector(dim) += projData(mode, count) * mainModeVects(dim, mode);
            }
          }

          index = 0;
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp = pOut + a * xoff + b * yoff + c * zoff;
                *pTemp += currVector(index);
                pTemp2 = pCounts + a * xoff + b * yoff + c * zoff;
                *pTemp2 += 1;
                ++index;
              }
            }
          }

          ++count;
        } // End of valid local windows section.

        ++pOut;
        ++pCounts;
        ++pMask;
      }
    }
  }

  //renormalise output image.
  pOut  = output->GetPointerToVoxels();
  pCounts = visitCounts->GetPointerToVoxels();
  for (k = radius; k < zmax; ++k){
    for (j = radius; j < ymax; ++j){
      for (i = radius; i < xmax; ++i){
        if (*pCounts > 0){
          *pOut = *pOut / ((double) *pCounts);
        }
        ++pOut;
        ++pCounts;
      }
    }
  }

  cout << "Data re-drawn using main mode basis." << endl;
}

void assignDists (irtkVector &dists, irtkGreyImage *mask, irtkRealImage *output, int radius){
  int i, j, k, count = 0;
  irtkGreyPixel *pMask;
  irtkRealPixel *pOut;

  // Limits.
  int xmax = mask->GetX() - 1 - radius;
  int ymax = mask->GetY() - 1 - radius;
  int zmax = mask->GetZ() - 1 - radius;
  // Main loop.
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pMask  = mask->GetPointerToVoxels(radius, j, k);
      pOut = output->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
        if(*pMask > 0){
          *pOut = dists(count);
          ++count;
        } else {
          *pOut = 0.0;
        }
        ++pMask;
        ++pOut;
      }

    }
  }
  cout << "distances assigned to image" << endl;
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
  cout << "normalised projected data" << endl;
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
  cout << "normalised distances calculated" << endl;
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
  irtkMatrix mainModeVects, projData, covMatrix;
  irtkVector meanWindow, eigenvals, dists;

  int radius = 1, i, j, ok;
  irtkGreyImage *input, *mask;
  irtkGreyPixel padding = MIN_GREY;
  irtkRealImage *output;

  int   windowCount, dims;
  int   mainModeCount = 0;
  float minPercentage = 1.0;
  float maxDist;

  // Parse arguments.
  if (argc < 4){
    usage();
  }

  input_name = argv[1];
  argv++; argc--;
  output_name = argv[1];
  argv++; argc--;
  radius      = atoi(argv[1]);
  argv++; argc--;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;      argv++;
      padding = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  _stepSize = radius;


  // Read and pre-process images.
  input  = new irtkGreyImage(input_name);
  mask   = new irtkGreyImage(input_name);
  output = new irtkRealImage(input_name);


  preProcess(input, padding);
  getMask(input, radius, mask, windowCount);

  // dims is the number of voxels in a window.
  dims = (int) pow(2 * radius + 1, 3.0);

  meanWindow.Initialize(dims);
  meanWindow *= 0.0;
  getMeanVector(input, mask, radius, dims, meanWindow);
  covMatrix.Initialize(dims, dims);
  getCovMatrix(covMatrix, input, mask, radius, dims, meanWindow);
  // Compute the eigen vals and vecs for the covariance matrix.
  irtkEigenAnalysis ea1(dims);
  for (i = 0; i < dims; ++i){
    for(j = 0; j < dims; ++j){
      ea1.Matrix(i, j) = covMatrix(i, j);  
    }
  }
  ea1.DecrSortEigenStuff();

  // How many modes contribute at least the min percentage of the variance?
  mainModeCount = countMainModes(ea1, dims, minPercentage);
  // Copy the main modes
  mainModeVects.Initialize(dims, mainModeCount);
  mainModeVects *= 0.0;
  for (i = 0; i < dims ;          ++i){
    for (j = 0; j < mainModeCount; ++j){
      mainModeVects(i, j) = ea1.Eigenvector(i, j);
    }
  }
  // Copy the main eigenvalues.
  eigenvals.Initialize(mainModeCount);
  for (i = 0; i < mainModeCount ; i++){
    eigenvals(i) = ea1.Eigenvalue(i);
  }

  // The data is projected onto a basis consisting of the main mode
  // vectors. The size is mainModeCount x the number of data items.
  projData.Initialize(mainModeCount, windowCount);
  getProjectedData(projData, mainModeVects, input, mask, 
                   mainModeCount, dims, radius, meanWindow);

  redrawImage (projData, mainModeVects,
               output, mask, mainModeCount,
               dims, radius, meanWindow );


  output->Write(output_name);



//   // Normalise each dimension based on its variance (e-value).
//   normaliseProjData(projData, eigenvals, mainModeCount, windowCount);

//   dists.Initialize(windowCount);
//   getDistances(dists, projData, mainModeCount, windowCount);
//   // Free memory
//   projData.Initialize(1, 1);

//   maxDist = getMaxDist(dists, windowCount);
//   if (maxDist > 0){
//     dists = dists * (100.0 / maxDist);
//   }

//   assignDists(dists, mask, output, radius);

//   output->Write(output_name);

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

/////////////////////////////////////////////////////////////
// // Also generates the mask.
// int getWindowCount(irtkGreyImage *image1, irtkGreyImage *image2, 
//                    int radius, irtkGreyImage *mask){
//   int i, j, k;
//   int a, b, c;
//   int count = 0;
//   int minVal1, maxVal1, minVal2, maxVal2;

//   irtkGreyPixel *pCentre1, *pCentre2, *pTemp1, *pTemp2;
//   irtkGreyPixel *pMask;
//   irtkGreyPixel high1, low1;
//   irtkGreyPixel high2, low2;

//   int xmax = image1->GetX() - 1 - radius;
//   int ymax = image1->GetY() - 1 - radius;
//   int zmax = image1->GetZ() - 1 - radius;
//   int xoff = 1;
//   int yoff = image1->GetX();
//   int zoff = image1->GetX() * image1->GetY();

//   image1->GetMinMax(&low1, &high1);
//   image2->GetMinMax(&low2, &high2);

//   for (k = radius; k < zmax; ++k){
//     for (j = radius; j < ymax; ++j){

//       pCentre1 = image1->GetPointerToVoxels(radius, j, k);
//       pCentre2 = image2->GetPointerToVoxels(radius, j, k);
//       pMask    = mask->GetPointerToVoxels(radius, j, k);

//       for (i = radius; i < xmax; ++i){

//         minVal1 = high1;
//         maxVal1 = low1;
//         minVal2 = high2;
//         maxVal2 = low2;

//         for (c = -1 * radius; c <= radius; ++c){
//           for (b = -1 * radius; b <= radius; ++b){
//             for (a = -1 * radius; a <= radius; ++a){
//               pTemp1 = pCentre1 + a * xoff + b * yoff + c * zoff;
//               pTemp2 = pCentre2 + a * xoff + b * yoff + c * zoff;

//               if (minVal1 > *pTemp1)
//                 minVal1 = *pTemp1;
//               if (maxVal1 < *pTemp1)
//                 maxVal1 = *pTemp1;
//               if (minVal2 > *pTemp2)
//                 minVal2 = *pTemp2;
//               if (maxVal2 < *pTemp2)
//                 maxVal2 = *pTemp2;
//             }
//           }
//         }

//         if ((minVal1 > 0 && maxVal1 > minVal1) ||
//             (minVal2 > 0 && maxVal2 > minVal2) ){
//           ++count;
//           *pMask = 1;
//         } else {
//           *pMask = 0;
//         }

//         ++pCentre1;
//         ++pCentre2;
//         ++pMask;
//       }
//     }
//   }

//   return count;
// }
