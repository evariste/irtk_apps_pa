#include <irtkImage.h>
#include <irtkHistogram.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>
#include <nr.h>

char *input_name1 = NULL, *output_name = NULL;
char *input_name2 = NULL;
int _stepSize = 1;
bool blnQuiet = false;

void usage(){

  cout << "histogramCCA [imageA] [imageB] [radius]" << endl;
  cout << "-TpA: padding A. -TpB: padding B." << endl;
  //  cout << "-output: name for output image containing histogram, if required." << endl;
  cout << "-q: Quiet, reduced output, does not print all scree data." << endl;
  cout << "-overlap: make windows overlap, default is tiled." << endl;
  cout << "-modeCount: How many modes to use from the canonical correlates." << endl;
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

        if ((minVal1 > 0 && maxVal1 > minVal1 ) ||
            (minVal2 > 0 && maxVal2 > minVal2 ) ){
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
void getMeanVector(irtkGreyImage *image1,
                   irtkGreyImage *image2,
                   irtkGreyImage *mask,
                   int radius, int dims,
                   irtkVector &meanWindow){
  int i, j, k;
  int a, b, c;
  int index, count = 0;

  irtkGreyPixel *pCentre1, *pCentre2, *pTemp1, *pTemp2, *pMask;

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
  // Main loop.  
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre1 = image1->GetPointerToVoxels(radius, j, k);
      pCentre2 = image2->GetPointerToVoxels(radius, j, k);
      pMask    =   mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
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
        }
        pCentre1 += _stepSize;
        pCentre2 += _stepSize;
        pMask += _stepSize;
      }
    }
  }

//   cout << "\nmean window: "<< endl;
//   meanWindow.Print();

  meanWindow = meanWindow * (1.0 / count);
}
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// This routine avoids the overhead of storing the data for all windows in
// a single (huge) data matrix.
void getCovMatrix(irtkMatrix &covMatrix,
                  irtkGreyImage *image1,
                  irtkGreyImage *image2,
                  irtkGreyImage *mask,
                  int radius, int dims,
                  irtkVector &meanWindow){
  int i, j, k;
  int a, b, c, d1, d2;
  int count = 0;
  int index;
  irtkGreyPixel *pCentre1, *pCentre2, *pTemp1, *pTemp2, *pMask;
  irtkMatrix tempMat;
  irtkVector currVector;
  tempMat.Initialize(dims, dims);
  currVector.Initialize(dims);
  // Limits.
  int xmax = image1->GetX() - 1 - radius;
  int ymax = image1->GetY() - 1 - radius;
  int zmax = image1->GetZ() - 1 - radius;
  // Offsets.
  int xoff = 1;
  int yoff = image1->GetX();
  int zoff = image1->GetX() * image1->GetY();

  // Main loop.
  for (k = radius; k < zmax; k += _stepSize){
    for (j = radius; j < ymax; j += _stepSize){

      pCentre1 = image1->GetPointerToVoxels(radius, j, k);
      pCentre2 = image2->GetPointerToVoxels(radius, j, k);
      pMask    = mask->GetPointerToVoxels(radius, j, k);

      for (i = radius; i < xmax; i += _stepSize){
        if(*pMask > 0){
          //Valid local windows.
          currVector *= 0.0;
          index = 0;
          //Loop over window in first image.
          for (c = -1 * radius; c <= radius; ++c){
            for (b = -1 * radius; b <= radius; ++b){
              for (a = -1 * radius; a <= radius; ++a){
                pTemp1 = pCentre1 + a * xoff + b * yoff + c * zoff;
                currVector(index) = *pTemp1;
                ++index;
              }
            }
          }
          //Loop over window in second image.
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

          // add (x - mu)^T * (x - mu) to the accumulating cov matrix.
          for (d1 = 0; d1 < dims; ++d1){
            for (d2 = d1; d2 < dims; ++d2){
              tempMat(d1, d2) = tempMat(d2, d1) = currVector(d1) * currVector(d2);
            }
          }
          covMatrix += tempMat;
        } // End of valid local windows section.

        pCentre1 += _stepSize;
        pCentre2 += _stepSize;
        pMask += _stepSize;
      }
    }
  }
  // Normalise the cov matrix.
  covMatrix = covMatrix * (1.0 / (count - 1));
}
/////////////////////////////////////////////////////////////////
float getScreePlotData(irtkEigenAnalysis &ea, int noRequired, int dims){

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

  for (i = 0; i < noRequired; i++){

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


bool equalImageSizes(irtkGreyImage *first, irtkGreyImage *second){
  return ( first->GetX() == second->GetX() ) &&
         ( first->GetY() == second->GetY() ) &&
         ( first->GetZ() == second->GetZ() );
}


int main(int argc, char **argv)
{
  irtkMatrix covMatrix, Cxx, Cxy, Cyx, Cyy, Cderived;
  irtkVector meanWindow;
  irtkVector eigenvals;

  int radius = 1, i, j, ok;
  irtkGreyImage *input1, *input2, *mask;
  irtkGreyPixel padding1 = MIN_GREY, padding2 = MIN_GREY;

  int   windowCount, dims, singleWindowLength;
  int modesRequired = -1;
  int overlapping = False;
  float screePlotArea = 0.0;

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
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-TpA") == 0)){
      argc--;      argv++;
      padding1 = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-TpB") == 0)){
      argc--;      argv++;
      padding2 = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-overlap") == 0)){
      overlapping = True;
      argv++; argc--;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      blnQuiet = true;
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-modeCount") == 0)){
      argc--;      argv++;
      modesRequired = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Decide if the windows are overlapping or tiled (default).
  if (overlapping == True){
    _stepSize = 1;
  } else {
    _stepSize = 2 * radius + 1;
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
  cout << "Window count : " << windowCount << endl;

  // dims is the number of voxels in two concatenated windows.
  singleWindowLength = (int) pow(2 * radius + 1, 3.0);
  dims = 2 * singleWindowLength;

  meanWindow.Initialize(dims);
  meanWindow *= 0.0;
  // Get mean vector from concatenated windows in both images.
  getMeanVector(input1, input2, mask, radius, dims, meanWindow);

  covMatrix.Initialize(dims, dims);
  getCovMatrix(covMatrix, input1, input2, mask, radius, dims, meanWindow);

  // Extract Cxx, Cxy, Cyx and Cyy.
  Cxx.Initialize(singleWindowLength, singleWindowLength);
  Cxy.Initialize(singleWindowLength, singleWindowLength);
  Cyx.Initialize(singleWindowLength, singleWindowLength);
  Cyy.Initialize(singleWindowLength, singleWindowLength);
  Cderived.Initialize(singleWindowLength, singleWindowLength);

  for (i = 0; i < singleWindowLength; ++i){
    for(j = 0; j < singleWindowLength; ++j){
      Cxx(i, j) = covMatrix(i         , j         );
      Cxy(i, j) = covMatrix(i         , j + singleWindowLength);
      Cyx(i, j) = covMatrix(i + singleWindowLength, j         );
      Cyy(i, j) = covMatrix(i + singleWindowLength, j + singleWindowLength);
    }
  }

  // The eigenvalues of the matrix Cderived form the canonical coefficients
  // of determination, Cderived = inv(Cxx)*Cxy*invCyy*Cyx .
  irtkMatrix invX, invY;
  invX.Initialize(singleWindowLength, singleWindowLength);
  invY.Initialize(singleWindowLength, singleWindowLength);

  invX = Cxx;  invX.Invert();
  invY = Cyy;  invY.Invert();

  Cderived = invX * Cxy * invY * Cyx;
  //Cderived = (!Cxx) * Cxy * (!Cyy) * Cyx;

  // Cderived is, generally, non-symmetric. So cannot use the
  // tridiagonalise and QL algorithm in irtkEigenAnalysis class. Use NR
  // routine appropriate for non-symmetric matrices instead.

  // Storage for NR routine.
  float **m = new float *[1 + singleWindowLength];
  for (i = 0; i < singleWindowLength + 1; ++i){
    m[i] = new float[singleWindowLength + 1];
  }

  // Eigenvalues, wr = real parts, wi = imaginary parts.
  float wr[singleWindowLength + 1];
  float wi[singleWindowLength + 1];

  // Numerical recipes arrays are 1-indexed.
  for (i = 0; i < singleWindowLength; ++i){
    for(j = 0; j < singleWindowLength; ++j){
      m[i+1][j+1] = Cderived(i, j);
    }
  }

  // Balancing stage.
  balanc(m , singleWindowLength);
  // Reduction to Hessenberg form (This form has the same eigenvalues as
  // the original matrix).
  elmhes(m,  singleWindowLength);
  // QR algorithm applied to Hessenberg form.
  hqr(m, singleWindowLength, wr, wi);

  sort2(singleWindowLength, wr, wi);

  if (!blnQuiet){
    for (i = 1; i <= singleWindowLength; ++i){
      cout << "i: " << i << "\t wr: " << wr[i] << "\t wi: " << wi[i] << endl;
    }
  }

// Checked the other way around ... i.e. used 
// Cderived = invY * Cyx * invX * Cxy; 
// and obtained same eigenvectors except for very small precision
// variations.

  float fTotal = 0.0, fCumulated = 0.0;
  float screeArea  = 0.0, percentage;
  for (i = 1; i <= singleWindowLength; ++i){
    fTotal +=  wr[i];
    if(wi[i] > 1E-05){
      cout << "Warning: Imaginary part of eigenvector exceeds 1E-05." << endl;
    }
  }
  // Eigenvectors are sorted in ascending order so read from the end.
  for (i = singleWindowLength; i >= 1; i--){
    percentage = 100.0 * wr[i] / fTotal;
    fCumulated += percentage;
    screeArea += fCumulated;
    if (!blnQuiet)
      cout << wr[i] << "\t" << percentage  << "\t"<< fCumulated << "\t"<< screeArea << endl;
  }

  cout << "Scree area: " << screeArea << "\t" <<screeArea / (100.0 * singleWindowLength)  << endl;

  for (i = 0; i < singleWindowLength; ++i){
    delete [] m[i];
  }
  delete [] m;
}

//   covMatrix.Write("cov.mat");
//   Cxx.Write("Cxx.mat");
//   Cyy.Write("Cyy.mat");
//   Cxy.Write("Cxy.mat");
//   Cyx.Write("Cyx.mat");
//   cout << "\nmean window: "<< endl;
//   meanWindow.Print();
 


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


// int countMainModes(irtkEigenAnalysis &ea, int dims, float minPercentage){

//   int i, mainModeCount = 0;
//   float fTotalVar   = 0.0;
//   //  float fCumulated = 0.0;

//   // Find sum of eigenvalues.
//   for (i = 0; i < dims; i++){
//     fTotalVar += ea.Eigenvalue(i);
//   }

//   // Identify the main modes.
//   for (i = 0; i < dims; i++){

//     if (100 * ea.Eigenvalue(i) / fTotalVar < minPercentage)
//       break;

//     mainModeCount++;
// //     cout<< "Mode:" << i+1 << " Eigenvalue:" << ea.Eigenvalue(i) << endl;
// //     fCumulated += 100 * ea.Eigenvalue(i) / fTotalVar;
// //     cout << " Mode [" << i+1 << "] explains " 
// //          << 100 * ea.Eigenvalue(i) / fTotalVar
// //          << " % (" 
// //          << fCumulated 
// //          << " %) of shape variance." << endl;
//   }

//   cout  << "mainModeCount :" << mainModeCount << endl;
//   return mainModeCount;
// }

// void getProjectedData (irtkMatrix &projData, irtkMatrix &mainModeVects,
//                        irtkGreyImage *image, int mainModeCount,
//                        int windowCount, int dims, int radius, irtkGreyImage *mask,
//                        irtkVector &meanWindow
// ){
//   int i, j, k;
//   int a, b, c, mode, dim;
//   int count = 0;
//   int index;
//   irtkGreyPixel *pCentre, *pTemp, *pMask;
//   irtkVector currVector;
//   currVector.Initialize(dims);

//   // Limits.
//   int xmax = image->GetX() - 1 - radius;
//   int ymax = image->GetY() - 1 - radius;
//   int zmax = image->GetZ() - 1 - radius;
//   // Offsets.
//   int xoff = 1;
//   int yoff = image->GetX();
//   int zoff = image->GetX() * image->GetY();

//   // Main loop.
//   for (k = radius; k < zmax; k += _stepSize){
//     for (j = radius; j < ymax; j += _stepSize){

//       pCentre = image->GetPointerToVoxels(radius, j, k);
//       pMask   = mask->GetPointerToVoxels(radius, j, k);

//       for (i = radius; i < xmax; i += _stepSize){
//         if(*pMask > 0){
//           //Valid local windows.
//           currVector *= 0.0;
//           index = 0;
//           for (c = -1 * radius; c <= radius; ++c){
//             for (b = -1 * radius; b <= radius; ++b){
//               for (a = -1 * radius; a <= radius; ++a){
//                 pTemp = pCentre + a * xoff + b * yoff + c * zoff;
//                 currVector(index) = *pTemp;
//                 ++index;
//               }
//             }
//           }

//           // Find x - mu
//           currVector -= meanWindow;

//           for (mode = 0; mode <  mainModeCount; ++mode){
//             projData(mode, count) = 0.0;
//             for (dim = 0; dim < dims; ++dim){
//               projData(mode, count) += mainModeVects(dim, mode) * currVector(dim);
//             }
//           }

//           ++count;
//         } // End of valid local windows section.

//         pCentre += _stepSize;
//         pMask += _stepSize;
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





// from old version of main


//   invX.Write("invX.mat");
//   invY.Write("invY.mat");
//   Cderived = (!Cxx);
//   Cderived = Cderived * Cxy;
//   Cderived = Cderived * (!Cyy);
//   Cderived = Cderived * Cyx;

  // Compute the eigen vals and vecs for Cderived.
//   irtkEigenAnalysis ea(singleWindowLength);
//   for (i = 0; i < singleWindowLength; ++i){
//     for(j = 0; j < singleWindowLength; ++j){
//       ea.Matrix(i, j) = Cderived(i, j);  
//     }
//   }
//   ea.DecrSortEigenStuff();

//   if (modesRequired == -1){
//     modesRequired = singleWindowLength;
//   }
  //  screePlotArea = getScreePlotData(ea, modesRequired ,singleWindowLength);

  //  cout << "scree plot area: " << screePlotArea << endl << endl;

  // Check the other way around ...
//   Cderived = (!Cyy);
//   Cderived = Cderived * Cyx;
//   Cderived = Cderived * (!Cxx);
//   Cderived = Cderived * Cxy;
 //  Cderived = invY * Cyx * invX * Cxy;

//   irtkEigenAnalysis ea2(singleWindowLength);
//   for (i = 0; i < singleWindowLength; ++i){
//     for(j = 0; j < singleWindowLength; ++j){
//       ea2.Matrix(i, j) = Cderived(i, j);
//     }
//   }
//   ea2.DecrSortEigenStuff();
 
  //  screePlotArea = getScreePlotData(ea2, modesRequired ,singleWindowLength);
  //  cout << "scree plot area: " << screePlotArea << endl;
