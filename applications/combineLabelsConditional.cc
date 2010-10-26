#include <irtkImage.h>

#define MAX_LABELS 67
#define ROOT_2_PI 2.50662827463100

char *query_input_name = NULL, *output_name = NULL, **input_names_anat = NULL, **input_names_seg = NULL;

void usage()
{
  cerr << "Usage: combineLabelConditional [queryLabel] [queryAnatImage] ";
  cerr << "[output] [N] [inputAnat1] .. [inputAnatN] [inputSeg1] ... [inputSegN]" << endl;
  cerr << "All input images are assumed to contain corresponding anatomies and labellings that" << endl;
  cerr << "are in register." << endl;
  cerr << "Given the intensities in [queryAnatImage] create an image where voxels that" << endl;
  cerr << "have the highest a posteriori probability of having label [queryLabel] are" << endl;
  cerr << "marked." << endl;
  cerr << "Optional:" << endl;
  cerr << "-kernel : standard deviation of gaussian used for Parzen estimation of class distribution."<< endl; 
  exit(1);
}

inline double gaussian(double mu, double sigma, double x)
{
  return exp( -1.0*(x - mu)*(x - mu) / (2.0*sigma*sigma)) / (sigma*ROOT_2_PI);
}


int lookupLabelIndex(short label, short labelsUsed[MAX_LABELS], int countLabelsUsed ){

  for (int i = 0; i < countLabelsUsed; ++i){
    if (label == labelsUsed[i]){
      return i;
    }
  }

  cout << "Error: Label " << label << " Not Found." << endl;
  exit(1);
  return -1;
}


int main(int argc, char **argv)
{
  short queryLabel, queryIntensity;
  int i, j, k, noInputImages, voxels, validVoxelCount;
  bool ok;

  irtkGreyImage query, mask, output, input;
  int *offsets;
  irtkGreyPixel *pIn, *pOut, *pMask, *pQuery;
  int countLabelsUsed;
  short labelsUsed[MAX_LABELS];
  bool bFound;
  short **labels;
  short **intensities;
  short **votes;
  double **lHoods;
  double **posteriors;
  // Width of gaussian used for Parzen windows.
  double sigma = 1.0;
  short *maxProbLabels;

  // Check command line
  if (argc < 7){
    usage();
  }

  queryLabel = atoi(argv[1]);
  argc--;
  argv++;
  query_input_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  noInputImages = atoi(argv[1]);
  argc--;
  argv++;

  input_names_anat = new char *[noInputImages];
  for (i = 0; i < noInputImages; i++){
    input_names_anat[i] = argv[1];
    argc--;
    argv++;
  }
  input_names_seg = new char *[noInputImages];
  for (i = 0; i < noInputImages; i++){
    input_names_seg[i] = argv[1];
    argc--;
    argv++;
  }

  // Parse any remaining arguments.
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-kernel") == 0)){
      argc--;
      argv++;
      sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read images
  query.Read(query_input_name);
  mask.Read(input_names_anat[0]);
  pMask = mask.GetPointerToVoxels();
  voxels = mask.GetNumberOfVoxels();
  memset ((void*)pMask, 0, voxels*sizeof(irtkGreyPixel));

  // Valid voxels are all voxels labeled with the query label by at least
  // one input image.  First, prepare mask:
  for (i = 0; i < noInputImages; ++i){
    input.Read(input_names_seg[i]);
    pIn  = input.GetPointerToVoxels();
    pMask = mask.GetPointerToVoxels();

    for (j = 0; j < voxels; ++j){
      if (*pIn == queryLabel){
        *pMask = 1;
      }
      ++pIn;
      ++pMask;
    }
  }
  // Now actually count ...
  validVoxelCount = 0;
  pMask = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pMask > 0){
      ++validVoxelCount;
    }
    ++pMask;
  }

//   cout << "validVoxelCount : " << validVoxelCount << endl;

  // ... and make note of offsets
  offsets = new int[validVoxelCount];
  j = 0;
  pMask = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pMask > 0){
      offsets[j] = i;
      ++j;
    }
    ++pMask;
  }

  // Read labels and intensities from valid voxels.
  // Each array is noInputImages x validVoxelCount
  labels      = new short*[noInputImages];
  intensities = new short*[noInputImages];
  for (i = 0; i < noInputImages; ++i){
    labels[i]      = new short[validVoxelCount];
    intensities[i] = new short[validVoxelCount];
  }
  // Segmented labels.
  for (i = 0; i < noInputImages; ++i){
    input.Read(input_names_seg[i]);
    pIn = input.GetPointerToVoxels();
    for (j = 0; j < validVoxelCount; ++j){
      labels[i][j] = *(pIn + offsets[j]);
    }
  }
  // Anatomic intensities.
  for (i = 0; i < noInputImages; ++i){
    input.Read(input_names_anat[i]);
    pIn = input.GetPointerToVoxels();
    for (j = 0; j < validVoxelCount; ++j){
      intensities[i][j] = *(pIn + offsets[j]);
    }
  }

  // Initialise array to store the labels used.
  // The query label occupies the first position in the array.
  countLabelsUsed = 1;
  labelsUsed[0] = queryLabel;
  // Blank out the remaining positions.
  for (i = 1; i < MAX_LABELS; ++i){
    labelsUsed[i] = -1;
  }
  // Loop over the valid voxels of each image.
  for (i = 0; i < noInputImages; ++i){
    for (j = 0; j < validVoxelCount; ++j){
      bFound = false;
      for (k = 0; k < countLabelsUsed; ++k){
        if (labels[i][j] == labelsUsed[k]){
          bFound = true;
          break;
        }
      }
      if (!bFound){
        labelsUsed[countLabelsUsed] = labels[i][j];
        ++countLabelsUsed;
      }
    }
  }

//   // A bit of communication.
//   cout << "countLabelsUsed : " << countLabelsUsed << endl;
//   for (i = 0; i < countLabelsUsed; ++i){
//     cout << "Label[" << i + 1 << "] = " << labelsUsed[i] << endl;
//     //    cout << " - index = " << lookupLabelIndex(labelsUsed[i], labelsUsed, countLabelsUsed)  << endl;
//   }

  // Initialise votes, likelihoods and posteriors.
  // These arrays are size validVoxelCount x countLabelsUsed
  votes = new short*[validVoxelCount];
  lHoods = new double*[validVoxelCount];
  posteriors = new double*[validVoxelCount];
  for (j = 0; j < validVoxelCount; ++j){
    votes[j] = new short[countLabelsUsed];
    memset((void*) votes[j], 0, countLabelsUsed*sizeof(short));
    lHoods[j] = new double[countLabelsUsed];
    memset((void*) lHoods[j], 0, countLabelsUsed*sizeof(float));
    posteriors[j] = new double[countLabelsUsed];
    memset((void*) posteriors[j], 0, countLabelsUsed*sizeof(float));
  }

  // Add up votes. votes[j][k] gives the votes that a class label, indexed
  // by k, receives at valid voxel j.  Dividing votes[j][k] by the number
  // of input images gives the prior probability of obtaining this class
  // label at the voxel.
  for (i = 0; i < noInputImages; ++i){
    for (j = 0; j < validVoxelCount; ++j){
      k = lookupLabelIndex(labels[i][j], labelsUsed, countLabelsUsed);
      ++votes[j][k];
    }
  }

//   // Debug
//   int unequivocalQueryVoxels = 0, rowSum;
//   for (j = 0; j < validVoxelCount; ++j){

//     rowSum = 0;
//     for (i = 0; i < countLabelsUsed; ++i){
//       if (votes[j][i] > (noInputImages / 2))
//         ++unequivocalQueryVoxels;
// //       rowSum += votes[j][i];
//     }
// //     if (rowSum != noInputImages){
// //       cout << "Row sum error in votes " << j << endl;
// //       exit(1);
// //     }
//   }
//   cout << "Voxels with unequivocal labelling: " << unequivocalQueryVoxels << endl;
//   // End Debug


  // Evaluate the likelihood at each valid voxel of getting the query
  // intensity given a particular label.
  // The distribution of intensities for each label is calculated by
  // centering a Gaussian window on each input image intensity where the
  // input image has this label.
  pQuery = query.GetPointerToVoxels();

  for (j = 0; j < validVoxelCount; ++j){
    queryIntensity = *(pQuery + offsets[j]);

    double sum = 0.0;
    for (k = 0; k < countLabelsUsed; ++k){
      // Does the current voxel have any votes for this label?
      if (votes[j][k] > 0){
        // Loop over images at this location.
        for (i = 0; i < noInputImages; ++i){
          // Did the current image vote for this label?
          if (labels[i][j] == labelsUsed[k]){
            // This image contributes a gaussian value to the likelihood.
            // The gaussian is centred on the intensity value for the input
            // image and is evaluated at the query intensity.
            lHoods[j][k] += gaussian(intensities[i][j], sigma, queryIntensity);
          }
        }
        // All contributions for this label at this voxel have been made.
        // Normalise by the number of contributing images.
        lHoods[j][k] /= votes[j][k];
      } else {
        // No images support the current label at this location.
        lHoods[j][k] = 0.0;
      }
      sum += lHoods[j][k];
    }
    for (k = 0; k < countLabelsUsed; ++k){
      lHoods[j][k] = 100 * lHoods[j][k] / sum;
    }
  }
  

  // Calculate the posteriors, the priors are directly proportional to the
  // votes so the votes are used instead.  Also division by a normalising
  // constant is ignored.
  for (j = 0; j < validVoxelCount; ++j){
    for (k = 0; k < countLabelsUsed; ++k){
      posteriors[j][k] = lHoods[j][k] * votes[j][k];
    }
  }

  // Identify the labels with maximum posterior probability.
  maxProbLabels = new short[validVoxelCount];
  float max;
  int maxProbLabelIndex;
  for (j = 0; j < validVoxelCount; ++j){
    max = -1000;
    maxProbLabelIndex = -1;
    for (k = 0; k < countLabelsUsed; ++k){
      if (max < posteriors[j][k]){
        max = posteriors[j][k];
        maxProbLabelIndex = k;
      }
    }
    maxProbLabels[j] = labelsUsed[maxProbLabelIndex];
  }

  // Write the output image. Assign 1 if the maximum probability label
  // matches the query label, 0 otherwise.
  output.Read(query_input_name);
  pOut = output.GetPointerToVoxels();
  voxels = output.GetNumberOfVoxels();
  memset((void*) pOut, 0, voxels * sizeof(irtkGreyPixel));

  for (j = 0; j < validVoxelCount; ++j){
    if (maxProbLabels[j] == queryLabel)
      *(pOut + offsets[j]) = 1;
  }

  output.Write(output_name);

  // Clean up.
  delete [] offsets;
  delete [] maxProbLabels;
  delete [] input_names_anat;
  delete [] input_names_seg;

  for (i = 0; i < noInputImages; ++i){
    delete [] labels[i];
    delete [] intensities[i];
  }
  delete [] labels;
  delete [] intensities;

  for (i = 0; i < validVoxelCount; ++i){
    delete [] votes[i];
    delete [] lHoods[i];
    delete [] posteriors[i];
  }
  delete [] votes;
  delete [] lHoods;
  delete [] posteriors;

}
