#include <irtkImage.h>
#include <irtkHistogram.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>

#define MAX_OBS 1000

char *input_name1 = NULL, *output_name = NULL;
char *input_name2 = NULL;
bool blnQuiet = false;

void usage(){

  cout << "screePlotArea [imageA] [imageB] [radius]" << endl;
  cout << "-TpA: padding A. -TpB: padding B." << endl;
  cout << "-q: Quiet, reduced output, does not print all scree data." << endl;
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
////////////////////////////////////////////////////////////////
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

  // Clear mask image.
  pMask  = mask->GetPointerToVoxels();
  int voxels = mask->GetNumberOfVoxels();
  for (i = 0; i < voxels; ++i){
    *pMask = 0;
    ++pMask;
  }
  //Limits.
  int xmax = image1->GetX() - 1 - radius;
  int ymax = image1->GetY() - 1 - radius;
  int zmax = image1->GetZ() - 1 - radius;
  //Offsets.
  int xoff = 1;
  int yoff = image1->GetX();
  int zoff = image1->GetX() * image1->GetY();

  image1->GetMinMax(&low1, &high1);
  image2->GetMinMax(&low2, &high2);

  for (k = radius; k < zmax; ++k){
    for (j = radius; j < ymax; ++j){

      pCentre1 = image1->GetPointerToVoxels(radius, j, k);
      pCentre2 = image2->GetPointerToVoxels(radius, j, k);
      pMask    = mask->GetPointerToVoxels(radius, j, k);

//       for (i = radius; i < xmax && count < MAX_OBS; ++i){
      for (i = radius; i < xmax; ++i){

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
        }
        ++pCentre1;
        ++pCentre2;
        ++pMask;
      }
    }
  }
  return count;
}
//////////////////////////////////////////////////////////////////
void getMeanVector(irtkGreyImage *image1, irtkGreyImage *image2,
                   irtkGreyImage *mask,
                   int radius, int dims,
                   irtkVector &meanWindow){
  int i, j, k;
  int a, b, c;
  int index, count = 0;

  irtkGreyPixel *pCentre1, *pTemp1, *pCentre2, *pTemp2, *pMask;

  int xmax = image1->GetX() - 1 - radius;
  int ymax = image1->GetY() - 1 - radius;
  int zmax = image1->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = image1->GetX();
  int zoff = image1->GetX() * image1->GetY();

  // Initialise mean.
  for (i = 0; i < dims; ++i)
    meanWindow(i) = 0.0;
  

  for (k = radius; k < zmax; ++k){
    for (j = radius; j < ymax; ++j){

      pCentre1 = image1->GetPointerToVoxels(radius, j, k);
      pCentre2 = image2->GetPointerToVoxels(radius, j, k);
      pMask    = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; ++i){
        if(*pMask > 0){
          index = 0;
          //Loop over window in first image.
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp1 = pCentre1 + a * xoff + b * yoff + c * zoff;
                meanWindow(index) += *pTemp1;
                ++index;
              }
            }
          }
          //Loop over window in second image.
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp2 = pCentre2 + a * xoff + b * yoff + c * zoff;
                meanWindow(index) += *pTemp2;
                ++index;
              }
            }
          }
          count++;
        } // if *pMask > 0
        ++pCentre1;
        ++pCentre2;
        ++pMask;
      }
    }
  }
  meanWindow = meanWindow * (1.0 / count);
}
//////////////////////////////////////////////////////////////////////////
// This routine avoids the overhead of storing the data for all windows in
// a single (huge) data matrix.
void getCovMatrix(irtkMatrix &covMatrix,
                  irtkGreyImage *image1, irtkGreyImage *image2,
                  irtkGreyImage *mask,
                  int radius, int dims,
                  irtkVector &meanWindow){
  int i, j, k;
  int a, b, c, d1, d2;
  int count = 0;
  int index;
  irtkGreyPixel *pCentre1, *pTemp1, *pCentre2, *pTemp2, *pMask;
  irtkMatrix temp;
  irtkVector currVector;
  temp.Initialize(dims, dims);
  currVector.Initialize(dims);
  // Limits.
  int xmax = image1->GetX() - 1 - radius;
  int ymax = image1->GetY() - 1 - radius;
  int zmax = image1->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = image1->GetX();
  int zoff = image1->GetX() * image1->GetY();

  for (k = radius; k < zmax; ++k){
    for (j = radius; j < ymax; ++j){

      pCentre1 = image1->GetPointerToVoxels(radius, j, k);
      pCentre2 = image2->GetPointerToVoxels(radius, j, k);
      pMask    = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; ++i){

        if(*pMask > 0){
          //Valid local windows.
          currVector *= 0.0;
          index = 0;
          // Loop over window in first image.
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp1 = pCentre1 + a * xoff + b * yoff + c * zoff;
                currVector(index) = *pTemp1;
                ++index;
              }
            }
          }
          // Loop over window in second image.
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp2 = pCentre2 + a * xoff + b * yoff + c * zoff;
                currVector(index) = *pTemp2;
                ++index;
              }
            }
          }
          ++count;
          // Find x - mu
          currVector -= meanWindow;
          
          // add (x - mu)T * (x - mu) to the accumulating cov matrix.
          for (d1 = 0; d1 < dims; ++d1){
            for (d2 = d1; d2 < dims; ++d2){
              temp(d1, d2) = temp(d2, d1) = currVector(d1) * currVector(d2);
            }
          }
          covMatrix += temp;
        } // End of valid local windows section.

        ++pCentre1;
        ++pCentre2;
        ++pMask;
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
/////////////////////////////////////////////////////////////////
float getScreePlotData(irtkEigenAnalysis &ea, int dims){

  int i;
  float fTotalVar   = 0.0;
  float fCumulated = 0.0;
  float percentage;
  float screeArea = 0.0;

  // Find sum of eigenvalues.
  for (i = 0; i < dims; i++){
    fTotalVar += ea.Eigenvalue(i);
  }

  if (!blnQuiet){
    cout<< "Mode,Eigenvalue,Explains,Cumulative" << endl;
  }

  for (i = 0; i < dims; i++){

    percentage  = 100 * ea.Eigenvalue(i) / fTotalVar;
    fCumulated += percentage;
    screeArea  += fCumulated;

    if (!blnQuiet){
      cout << i+1              << ",";
      cout << ea.Eigenvalue(i) << ",";
      cout << percentage       << ",";
      cout << fCumulated       << endl;
    }
  }

  return screeArea / 100.0;
}
////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

  irtkMatrix dataMatrix, covMatrix;
  irtkVector meanWindow;
  irtkVector eigenvals;
  int windowCount, dims;
  int radius = 1, i, j;
  bool ok;
  irtkGreyImage *input1, *input2, *mask;
  irtkGreyPixel padding1 = MIN_GREY, padding2 = MIN_GREY;

  // Parse arguments.
  if (argc < 4){
    usage();
  }

  input_name1 = argv[1];
  argv++; argc--;
  input_name2 = argv[1];
  argv++; argc--;
  radius = atoi(argv[1]);
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
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      blnQuiet = true;
      argc--;      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
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

  // Corresponding windows are concatenated into a single vector.
  // dims is therefore number of voxels in two windows.
  dims = 2 * (int) pow(2 * radius + 1, 3.0);

  meanWindow.Initialize(dims);
  meanWindow *= 0.0;
  getMeanVector(input1, input2, mask, radius, dims, meanWindow);

  covMatrix.Initialize(dims, dims);
  getCovMatrix(covMatrix, input1, input2, mask, radius, dims, meanWindow);

  irtkEigenAnalysis ea(dims);
  for (i = 0; i < dims; ++i){
    for(j = 0; j < dims; ++j){
      ea.Matrix(i, j) = covMatrix(i, j);
    }
  }
  ea.DecrSortEigenStuff();

  cout << "Scree plot data: " << endl;
  float screeArea = getScreePlotData(ea, dims);
  cout << "Scree Plot Area: " << screeArea;
  cout << " (" << dims << " dimensions)" << endl;
  cout << "Normalised: " << screeArea / dims << endl;

}

//////////////////////////////////////////////////////////////////////////////////////

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


// void getProjectedData (irtkMatrix &projData, irtkMatrix &mainModeVects,
//                       irtkMatrix &dataMatrix, int mainModeCount, int windowCount, int dims){
//   int i, j, k;

//   for (i = 0; i < mainModeCount; ++i){
//     for (j = 0; j < windowCount; ++j){
//       projData(i, j) = 0.0;
//       for (k = 0; k < dims; ++k){
//         projData(i, j) += mainModeVects(k, i) * dataMatrix(k, j);
//       }
//     }
//   }
// }

// void normaliseProjData (irtkMatrix &projData, irtkVector &eigenvals, 
//                         int mainModeCount, int dataItems){
//   // projData contains the components of the original data w.r.t. the main
//   // eigen modes.
//   // Divide each component by the corresponding eigenvalue (variance).
//   // This is done to get the Mahalanobis distance later.
//   int i, j;
//   float variance;

//   for (i = 0; i < mainModeCount; ++i){
//     variance = eigenvals(i);
//     for (j = 0; j < dataItems;     ++j){
//       projData(i, j) = projData(i, j) / variance;
//     }
//   }
// }

// void getDistances (irtkVector &dists, irtkMatrix &projData, int mainModeCount, int dataItems){
//   int i, j;

//   for (j = 0; j < dataItems;     ++j){
//     dists(j) = 0.0;
//     for (i = 0; i < mainModeCount; ++i){
//       dists(j) += projData(i, j) * projData(i, j);
//     }
//     dists(j) = sqrt(dists(j));
//   }
// }

// float getMaxDist(irtkVector &dists, int dataItems){
//   int i;
//   float max = FLT_MIN;
//   for (i = 0; i < dataItems; ++i){
//     if (max < dists(i)){
//       max = dists(i);
//     }
//   }
//   return max;
// }
