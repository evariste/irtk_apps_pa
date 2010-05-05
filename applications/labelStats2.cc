#include <irtkImage.h>

#include <map.h>

char *imageA_name = NULL;
char *imageB_name = NULL;
char *labelFile   = NULL;

#define MAX_LABELS 1000


map<short, short>::iterator iter;


void usage()
{
  cout << "labelStats [imageA] [imageB] <-options>" << endl;
  cout << "Images are label images (max 1000)." << endl;
  cout << "Print, for each label, the number of voxels assigned"  << endl;
  cout << "that label by each separate image and the number of "  << endl;
  cout << "voxels assigned the label by both images."  << endl;
  cout << "i.e. n(A), n(B) and n(A&B)"  << endl;
  cout << "Images must have same dimensions."  << endl;
  cout << "Output format (comma separated):" << endl;
  cout << "Label,n(A),n(B),n(A^B),OR,SI" << endl;
  cout << "options:" << endl;
  cout << "-summary : print only the mean SI and OR values." << endl;
  cout << "-file    : File containg a list of labels over which the " << endl;
  cout << "           attention is restricted (separated by whitespace)." << endl;
  cout << "-false   : print false positives / negatives (verbose output only)" << endl;
  cout << "           format becomes:" << endl;
  cout << "           Label,n(A),n(B),n(A^B),n(A\\B),n(B\\A),OR,SI" << endl;
  cout << "-q       : quiet summary, csv output : 'OR,SI' " << endl;
  cout << "-siRow   : Comma separated SI values in a row only (no labels)." << endl;
  cout << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, nVoxels;
  irtkGreyPixel *pA, *pB, max, min, tmp;
  int jointLabels[MAX_LABELS][MAX_LABELS];
  int marginalA[MAX_LABELS];
  int marginalB[MAX_LABELS];
  int includeLabel[MAX_LABELS];

  int printSummary = False;
  int quiet = False;
  int printFalseValues = False;
  int ok, temp;
  int siRow = False;

  // Read arguments
  if (argc < 3){
    usage();
  }

  imageA_name = argv[1];
  argc--;
  argv++;
  imageB_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-summary") == 0)){
      argc--;
      argv++;
      printSummary = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-siRow") == 0)){
      argc--;
      argv++;
      siRow = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-file") == 0)){
      argc--;
      argv++;
      labelFile = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-false") == 0)){
      argc--;
      argv++;
      printFalseValues = True;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  if (quiet == True){
    // Force summary.
    printSummary = True;
  }


  irtkGreyImage *imageA, *imageB;
  imageA = new irtkGreyImage(imageA_name);
  imageB = new irtkGreyImage(imageB_name);

  nVoxels = imageA->GetNumberOfVoxels();

  if (imageA->GetX() != imageB->GetX() ||
      imageA->GetY() != imageB->GetY() || 
      imageA->GetZ() != imageB->GetZ()){
    cout << "Equal Sizes Please!" << endl;
    exit(1);
  }

  // Converting labels to indices:
  map <short, short> lab2ind;
  int currIndex;

  pA = imageA->GetPointerToVoxels();

  currIndex = 1;

  for (i = 0; i < nVoxels; ++i){
    if (*pA < 1){
      ++pA;
      continue;
    }

    if (lab2ind.count(*pA) < 1){
      lab2ind[*pA] = currIndex;
      ++currIndex;
    }

    ++pA;
  }

//   if (!quiet)
//     cout << "Found " << currIndex - 1 << " distinct labels in image A." << endl;

  temp = currIndex;

  pB = imageB->GetPointerToVoxels();
  for (i = 0; i < nVoxels; ++i){
    if (*pB < 1){
      ++pB;
      continue;
    }

    if (lab2ind.count(*pB) < 1){
      lab2ind[*pB] = currIndex;
      ++currIndex;
    }

    ++pB;
  }

//   if (!quiet)
//     cout << "Found " << currIndex - temp << " further distinct labels in image B." << endl;

  max = currIndex;

  if (max >= MAX_LABELS){
    cout << "Max number of distinct labels exceeds limit" << endl;
    exit(1);
  }

  // Reverse lookup, indices back to labels:
  map <short, short> ind2lab;
  for (iter = lab2ind.begin(); iter != lab2ind.end(); ++iter){
    ind2lab[iter->second] = iter->first;
  }

  //Reset histograms.
  for (i = 0; i < MAX_LABELS; ++i){
    marginalA[i] = 0;
    marginalB[i] = 0;
    for (j = 0; j < MAX_LABELS; ++j){
      jointLabels[i][j] = 0;
    }
  }

  // Reset the list of included labels.
  if (labelFile == NULL){
    for (i = 0; i < MAX_LABELS; ++i){
      includeLabel[i] = 1;
    } 
  } else {
    // First clear all the include labels.
    for (i = 0; i < MAX_LABELS; ++i){
      includeLabel[i] = 0;
    }

    // Read in those to be included.
    ifstream fileIn;
    fileIn.open(labelFile, ifstream::in);
    fileIn.flags(ios_base::skipws);
    while (!fileIn.eof()){
      fileIn >> temp;
      temp = lab2ind[temp];
      if (temp > 0 && temp < MAX_LABELS){
        includeLabel[temp] = 1;
      }
    }
    fileIn.close();
  }


  pA = imageA->GetPointerToVoxels();
  pB = imageB->GetPointerToVoxels();

  int indA, indB;

  // Fill in histogram.
  for (i = 0; i < nVoxels; ++i){
    if (*pA > 0 || *pB > 0){
      indA = lab2ind[*pA];
      indB = lab2ind[*pB];
      ++marginalA[indA];
      ++marginalB[indB];
      ++jointLabels[indA][indB];
    }
    ++pA;
    ++pB;
  }

  double overlap, si;
  double sumOverlap, sumSi;
  int count;

  if (printSummary == True){
    count = 0;
    sumOverlap = sumSi = 0;

    for (i = 1; i <= max; ++i){
      if (includeLabel[i] && (marginalA[i] + marginalB[i] > 0)){
        overlap = jointLabels[i][i] / ((double) ( marginalA[i] + marginalB[i] - jointLabels[i][i] ));
        si = 2 *  jointLabels[i][i] / ((double) marginalA[i] + marginalB[i] );
        sumOverlap += overlap;
        sumSi += si;
        ++count;
      }
    }

    if (quiet){
      if (count > 0){
        cout << sumOverlap / ((double) count) << ",";
        cout << sumSi / ((double) count) << endl;
      } else {
        cout << "0,0" << endl;
      }
    } else {
      if (count > 0){
        cout << "Mean OR : " << sumOverlap / ((double) count) << endl;
        cout << "Mean SI : " << sumSi / ((double) count) << endl;
      } else {
        cout << "Zero overlap between images." << endl;
      }
    }

  } else {

    // All labels.
    for (iter = lab2ind.begin(); iter != lab2ind.end(); ++iter){

//       cout << iter->first << " " << iter->second << endl;

      i = iter->second;

      if (iter->first < 1)
        continue;

//     for (i = 1; i <= max; ++i){
      if (includeLabel[i] && (marginalA[i] + marginalB[i] > 0)){

        overlap = jointLabels[i][i] / ((double) ( marginalA[i] + marginalB[i] - jointLabels[i][i] ));
        si = 2 *  jointLabels[i][i] / ((double) marginalA[i] + marginalB[i] );

        if (siRow == True){
          cout << si << ",";
        } else {
          cout << iter->first << "," << marginalA[i] << "," << marginalB[i] << "," << jointLabels[i][i];
          if (printFalseValues == True){
            cout << "," << marginalA[i] - jointLabels[i][i] << "," << marginalB[i] - jointLabels[i][i];
          }
          cout << "," << overlap << "," << si << endl;
        }

      }
    }

    if (siRow == True){
      cout << endl;
    }
  }

}
