#include <irtkImage.h>

char *anat_name  = NULL,
  *probsA_name = NULL, *probsB_name = NULL, *probsC_name = NULL,
  *output_name = NULL;

#define PROB_RANGE 256

void usage();

void generateMask(irtkGreyImage &probsA, irtkGreyImage &probsB,
                  irtkGreyImage &probsC, irtkGreyImage &mask,
                  int total);

void getParams(double &mu, double &var, double &vol, int total,
               irtkGreyImage *probsMap, irtkGreyImage *anat, irtkGreyImage *mask);

int main(int argc, char **argv)
{

  irtkGreyImage anat;
  irtkGreyPixel *pA, *pB, *pC, *pMask, *pOut;
  irtkGreyImage probsA, probsB, probsC;
  irtkGreyImage output;
  irtkGreyImage probsOther;
  irtkGreyImage mask;
  //  irtkGreyPixel  minInt, maxInt;
  int i, voxels; // , intRange, index;
  double sum;
  irtkGreyPixel tol, total;

  //cluster parameters
  double muA,  muB,  muC; //,  muOther;
  double varA, varB, varC; //, varOther;
  double volA, volB, volC; //, volOther, volTotal;
  //  double wgtA, wgtB, wgtC;
  int meanA, meanB, meanC;

  if (argc < 6){
    usage();
  }

  anat_name   = argv[1];  argc--;  argv++;
  probsA_name = argv[1];  argc--;  argv++;
  probsB_name = argv[1];  argc--;  argv++;
  probsC_name = argv[1];  argc--;  argv++;
  output_name = argv[1];  argc--;  argv++;

  total = atoi(argv[1]);  argc--;  argv++;
  tol   = atoi(argv[1]);

  anat.Read(anat_name);
  output.Read(anat_name);
  mask.Read(anat_name);

  probsA.Read(probsA_name);
  probsB.Read(probsB_name);
  probsC.Read(probsC_name);

  generateMask(probsA, probsB, probsC, mask, total);

  getParams(muA, varA, volA, total, &probsA, &anat, &mask);
  getParams(muB, varB, volB, total, &probsB, &anat, &mask);
  getParams(muC, varC, volC, total, &probsC, &anat, &mask);

  meanA = (int) round(muA);
  meanB = (int) round(muB);
  meanC = (int) round(muC);

  voxels = probsA.GetNumberOfVoxels();

  pMask = mask.GetPointerToVoxels();
  pA    = probsA.GetPointerToVoxels();
  pB    = probsB.GetPointerToVoxels();
  pC    = probsC.GetPointerToVoxels();
  pOut  = output.GetPointerToVoxels();

  cout << "Performing estimates ... " << endl;

  for (i = 0; i < voxels; ++i){
    if (*pMask > 0){
      if (*pC < tol){
        if (*pB < tol){
          // 1 wm voxel 
          *pOut = meanA;
        } else if (*pA < tol){
          // 2 gm voxel
          *pOut = meanB;
        } else { 
          // 3 mixed wm gm
          sum = *pA + *pB;
          *pOut = (int) round( ( (*pA) * meanA + (*pB) * meanB ) / sum  );
        }
      } else { 
        //csf significant
        if (*pA < tol){
          if (*pB < tol){
            // 4 csf voxel
            *pOut = meanC;
          } else {
            // 5 gm csf mix
            sum = *pB + *pC;
            *pOut = (int) round( ( (*pB) * meanB + (*pC) * meanC ) / sum  );
          }
        } else {
          // wm significant
          if (*pB < tol){
            // 6 wm csf mix??
            sum = *pA + *pC;
            *pOut = (int) round( ( (*pA) * meanA + (*pC) * meanC ) / sum  );
          } else {
            // wm gm csf mix??
            sum = *pA + *pB + *pC;
            *pOut = (int) round( ( (*pA) * meanA + (*pB) * meanB + (*pC) * meanC ) / sum  );
          }
        }
      }
    }else {
      *pOut = 0;
    }


    pMask++;
    pA++;
    pB++;
    pC++;
    pOut++;
  }

  output.Write(output_name);

}

void getParams(double &mu, double &var, double &vol, int total,
               irtkGreyImage *probsMap, irtkGreyImage *anat, irtkGreyImage *mask){

  int i, voxels;
  irtkGreyPixel *pProb;
  irtkGreyPixel *pAnatInt, *pMask;
  double sum = 0.0, sumSquares = 0.0;
  double temp;

  cout << "Estimating cluster parameters ..." << endl;

  voxels = probsMap->GetNumberOfVoxels();

  pProb    = probsMap->GetPointerToVoxels();
  pAnatInt = anat->GetPointerToVoxels();
  pMask    = mask->GetPointerToVoxels();

  vol = mu = var = 0.0;

  for (i = 0; i < voxels; ++i){

    if (*pMask > 0){
      temp = ((double) *pProb) / total;
      vol        += temp;
      sum        += temp * (*pAnatInt);
      sumSquares += temp * (*pAnatInt) * (*pAnatInt);
    }
    ++pProb;
    ++pAnatInt;
    ++pMask;
  }

  mu  = sum / vol;
  var = (sumSquares / vol) - (mu * mu);

  cout << "done" << endl;
  cout << "Mean   = " << mu << endl;
  cout << "Var    = " << var << endl;
  cout << "Volume = " << vol << endl;
}

// Also finds the probability of class 'other'.
void generateMask(irtkGreyImage &probsA, irtkGreyImage &probsB,
                  irtkGreyImage &probsC, irtkGreyImage &mask,
                  int total){
  irtkGreyPixel *pA, *pB, *pC, *pMask;
  int i, voxels;

  cout << "Generating mask ... " << endl;

  // Pre: all images same size.
  voxels = probsA.GetNumberOfVoxels();
  pA     = probsA.GetPointerToVoxels();
  pB     = probsB.GetPointerToVoxels();
  pC     = probsC.GetPointerToVoxels();
  pMask  = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pA + *pB + *pC < total){
      *pMask  = 0;
    } else {
      *pMask  = 100;
    }

    ++pA;    ++pB;    ++pC;
    ++pMask;
  }
  cout << "done" << endl;
}



void usage()
{
  cerr << "\t Usage:" << endl;
  cerr << "\t simulateMRimage [anatImage] [probsWM] [probsGM] [probsCSF] [output] [total] [tol]" << endl;
  cerr << "\t Generate a simulated MR image from a real anatomical image and" << endl;
  cerr << "\t and three images representing the probability maps for the tissue" << endl;
  cerr << "\t tissue classes (WM, GM, CSF)." << endl;
  cerr << "\t Any bias field is not taken into account." << endl;
  cerr << "\t [Total] is what the probabilities add up to, [tol] is the level to " << endl;
  cerr << "\t decide that a particular tissue is *not* present in a voxel (i.e. lower limit)." << endl;

  exit(1);
}

/*

////////////////////////////////////////////////////////////////////////
/// VERSION 1

#include <irtkImage.h>

#include "vtkIdList.h"
#include "vtkPoints.h"
#include <vtkPointLocator.h>
#include <vtkUnstructuredGrid.h>

char *anat_name  = NULL, 
     *probsA_name = NULL, *probsB_name = NULL, *probsC_name = NULL, 
     *output_name = NULL;

#define PROB_RANGE 256

void usage();

void normalise(irtkGreyImage &img);

// Also finds the probability of class 'other'.
void generateMask(irtkGreyImage &probsA, irtkGreyImage &probsB,
                  irtkGreyImage &probsC, irtkGreyImage &probsOther,
                  irtkGreyImage &mask);

void getParams(double &mu, double &var, double &vol, 
               irtkGreyImage *probsMap, irtkGreyImage *anat, irtkGreyImage *mask);

inline double gaussian(double mu, double var, double x){

  return exp(-1.0 * (x - mu) * (x - mu) / (2.0 * var)) / sqrt(2.0 * M_PI * var);

}

void calcLookup(double *likelihoods, double mu, double var, 
                int range, int minInt, double weight);


int main(int argc, char **argv)
{

  irtkGreyImage anat;
  irtkGreyImage probsA, probsB, probsC;
  irtkGreyImage output;
  irtkGreyImage probsOther;
  irtkGreyImage mask;
  irtkGreyPixel  minInt, maxInt;  
  int i, tol, voxels, intRange, index;
  double sum;

  // pixel pointers
  irtkGreyPixel *pMask, *pProbsA, *pProbsB, *pProbsC, *pOut;

  //cluster parameters
  double muA,  muB,  muC,  muOther;
  double varA, varB, varC, varOther;
  double volA, volB, volC, volOther, volTotal;
  double wgtA, wgtB, wgtC;
  double *likelihoodsA, *likelihoodsB,*likelihoodsC;
  double posteriors[3];

  if (argc < 6){
    usage();
  }

  anat_name   = argv[1];  argc--;  argv++;
  probsA_name = argv[1];  argc--;  argv++;
  probsB_name = argv[1];  argc--;  argv++;
  probsC_name = argv[1];  argc--;  argv++;
  output_name = argv[1];  argc--;  argv++;
  tol         = atoi(argv[1]);

  anat.Read(anat_name);
  output.Read(anat_name);

  probsA.Read(probsA_name);
  probsB.Read(probsB_name);
  probsC.Read(probsC_name);

  // These images will be over-written during mask generation.
  probsOther.Read(probsA_name);
  mask.Read(anat_name);

  normalise(probsA);
  normalise(probsB);
  normalise(probsC);

  generateMask(probsA, probsB, probsC, probsOther, mask);

  getParams(muA, varA, volA, &probsA, &anat, &mask);
  getParams(muB, varB, volB, &probsB, &anat, &mask);
  getParams(muC, varC, volC, &probsC, &anat, &mask);
  getParams(muOther, varOther, volOther, &probsOther, &anat, &mask);

//   probsA.Write("probsA.hdr");
//   probsB.Write("probsB.hdr");
//   probsC.Write("probsC.hdr");
//   probsOther.Write("probsOther.hdr");

  volTotal = volA + volB + volC + volOther;
  if (volOther / volTotal > 0.005){
    cout << "Volume 'other' is " << 100 * volOther / volTotal << "% of the total" << endl;
  }

  //weights
  wgtA = volA / volTotal;
  wgtB = volB / volTotal;
  wgtC = volC / volTotal;
  cout << "Weights: " << wgtA << "\t"  << wgtB << "\t"  << wgtC << endl;

  anat.GetMinMax(&minInt, &maxInt);
  cout << "Min and max intensities: " << minInt << "\t" << maxInt << endl;
  intRange = maxInt - minInt + 1;

  likelihoodsA = new double[intRange];
  likelihoodsB = new double[intRange];
  likelihoodsC = new double[intRange];

  calcLookup(likelihoodsA, muA, varA, intRange, minInt, wgtA);
  calcLookup(likelihoodsB, muB, varB, intRange, minInt, wgtB);
  calcLookup(likelihoodsC, muC, varC, intRange, minInt, wgtC);

  // populate the point locator, each point is a triple of probabilities:
  vtkPoints     *points  = vtkPoints::New();

  cout << "Setting up locator ..." << endl;

  for (i = 0; i < intRange; ++i){
    posteriors[0] = likelihoodsA[i];
    posteriors[1] = likelihoodsB[i];
    posteriors[2] = likelihoodsC[i];

    sum = posteriors[0] + posteriors[1] + posteriors[2];

    posteriors[0] = posteriors[0] / sum;
    posteriors[1] = posteriors[1] / sum;
    posteriors[2] = posteriors[2] / sum;

    //    cout << i << "\t" << posteriors[0] << "\t" << posteriors[1] << "\t" << posteriors[2] << endl;
    points->InsertPoint(i, &posteriors[0]);
  }

  vtkUnstructuredGrid *grid = vtkUnstructuredGrid::New();
  grid->SetPoints(points);

  vtkPointLocator *pLoc = vtkPointLocator::New();
  pLoc->SetDataSet(grid);
  pLoc->BuildLocator();

  // Find closest match based on probability maps given.
  voxels = probsA.GetNumberOfVoxels();

  pMask   = mask.GetPointerToVoxels();
  pProbsA = probsA.GetPointerToVoxels();
  pProbsB = probsB.GetPointerToVoxels();
  pProbsC = probsC.GetPointerToVoxels();
  pOut    = output.GetPointerToVoxels();

  cout << "Performing lookups ... " << endl;

  for (i = 0; i < voxels; ++i){
    sum = *pProbsA + *pProbsB + *pProbsC;

    if ((*pMask > 0) && 
        ( sum > tol)){
      //index into intensities, also indices of lookups
      index = pLoc->FindClosestPoint((double) *pProbsA / sum,
                                     (double) *pProbsB / sum,
                                     (double) *pProbsC / sum );
      *pOut = index + minInt;
    } else {

      *pOut = 0;
    }

    pMask++;
    pProbsA++;
    pProbsB++;
    pProbsC++;
    pOut++;

  }

  output.Write(output_name);
  //  i = pLoc->FindClosestPoint(a, b, c);
  //  points->GetPoint(i, pt);

  delete [] likelihoodsA;
  delete [] likelihoodsB;
  delete [] likelihoodsC;

}


void usage()
{
  cerr << "\t Usage:" << endl;
  cerr << "\t simulateMRimage [anatImage] [probsA] [probsB] [probsC] [output] [tol]" << endl;
  cerr << "\t Generate a simulated MR image from a real anatomical image and" << endl;
  cerr << "\t and three images representing the probability maps for the tissue" << endl;
  cerr << "\t tissue classes (WM, GM, CSF)." << endl;
  cerr << "\t Any bias field is not taken into account." << endl;
  cerr << "\t Images representing probabilities will be normalised to the range" << endl;
  cerr << "\t [0-256] so must not contain any values smaller than 0." << endl;
  cerr << "\t " << endl;
  cerr << "\t " << endl;

  exit(1);
}

void normalise(irtkGreyImage &img){

  // Assumes that the max value corresponds to a probability of 1,
  // i.e. that there are some voxels with an unequivocal probability set.
  cout << "normalising image ..." << endl;

  irtkGreyPixel  max = MIN_GREY;
  irtkGreyPixel *pPix;
  int i, voxels;

  voxels = img.GetNumberOfVoxels();
  pPix   = img.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pPix > max){
      max = *pPix;
    }
    ++pPix;
  }

  if (max <= 0){
    cerr << "Error normalising image, max <= 0." << endl;
    exit(1);
  }

  pPix = img.GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    *pPix = (int) round(PROB_RANGE * (*pPix / (double) max));
    ++pPix;
  }

  cout << "done" << endl;
}

// Also finds the probability of class 'other'.
void generateMask(irtkGreyImage &probsA, irtkGreyImage &probsB,
                  irtkGreyImage &probsC, irtkGreyImage &probsOther,
                  irtkGreyImage &mask){

  cout << "generate mask and 'other' probs ... " << endl;

  int tol = 1;
  irtkGreyPixel *pA, *pB, *pC, *pOther;
  irtkGreyPixel *pMask;
  int i, voxels;

  // Pre: all images same size.
  voxels = probsA.GetNumberOfVoxels();

  pA     = probsA.GetPointerToVoxels();
  pB     = probsB.GetPointerToVoxels();
  pC     = probsC.GetPointerToVoxels();
  pOther = probsOther.GetPointerToVoxels();
  pMask  = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){

    if (*pA < tol && *pB < tol && *pC < tol){
      *pMask  = 0;
      *pOther = 0;
    } else {
      *pMask  = 100;
      *pOther = PROB_RANGE - *pA - *pB - *pC;

      if (*pOther < 0){
        *pOther = 0;
      }
    }

    ++pA;    ++pB;    ++pC;
    ++pOther;
    ++pMask;

  }
  cout << "done" << endl;
}

void getParams(double &mu, double &var, double &vol, 
               irtkGreyImage *probsMap, irtkGreyImage *anat, irtkGreyImage *mask){

  int i, voxels;
  irtkGreyPixel *pProb;
  irtkGreyPixel *pAnatInt, *pMask;
  double sum = 0.0, sumSquares = 0.0;
  double temp;

  cout << "Estimating cluster parameters ..." << endl;

  voxels = probsMap->GetNumberOfVoxels();

  pProb    = probsMap->GetPointerToVoxels();
  pAnatInt = anat->GetPointerToVoxels();
  pMask    = mask->GetPointerToVoxels();

  vol = mu = var = 0.0;

  for (i = 0; i < voxels; ++i){

    if (*pMask > 0){
      temp = ((double) *pProb) / PROB_RANGE;
      vol        += temp;
      sum        += temp * (*pAnatInt);
      sumSquares += temp * (*pAnatInt) * (*pAnatInt);
    }
    ++pProb;
    ++pAnatInt;
    ++pMask;
  }

  mu  = sum / vol;
  var = (sumSquares / vol) - (mu * mu);

  cout << "done" << endl;
  cout << "Mean   = " << mu << endl;
  cout << "Var    = " << var << endl;
  cout << "Volume = " << vol << endl;
}

void calcLookup(double *likelihoods, double mu, double var, 
                int range, int minInt, double weight){

  for (int i = 0; i < range; ++i){
    likelihoods[i] = weight * gaussian(mu, var, i + minInt);
  }
}

*/
