#include <irtkImage.h>

char *imageA_name = NULL;
char *imageB_name = NULL;
char *labelFile   = NULL;

#define MAX_LABELS 1000

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
  int orRow = False;

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
    if ((ok == False) && (strcmp(argv[1], "-orRow") == 0)){
      argc--;
      argv++;
      orRow = True;
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

  if (imageA->GetX() != imageB->GetX() ||
      imageA->GetY() != imageB->GetY() || 
      imageA->GetZ() != imageB->GetZ()){
    cout << "Equal Sizes Please!" << endl;
    exit(1);
  }

  imageA->GetMinMax(&min, &max);
  tmp = max;
  imageB->GetMinMax(&min, &max);
  if (tmp > max)
    max = tmp;

  if (max >= MAX_LABELS){
    cout << "Max label exceeds limit" << endl;
    exit(1);
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
      if (temp > 0 && temp < MAX_LABELS){
        includeLabel[temp] = 1;
      }
    }
    fileIn.close();
  }


  pA = imageA->GetPointerToVoxels();
  pB = imageB->GetPointerToVoxels();
  nVoxels = imageA->GetNumberOfVoxels();

  // Fill in histogram.
  for (i = 0; i < nVoxels; ++i){
    if (*pA > 0 || *pB > 0){
      ++marginalA[*pA];
      ++marginalB[*pB];
      ++jointLabels[*pA][*pB];
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
    for (i = 1; i <= max; ++i){
      if (includeLabel[i] && (marginalA[i] + marginalB[i] > 0)){

        overlap = jointLabels[i][i] / ((double) ( marginalA[i] + marginalB[i] - jointLabels[i][i] ));
        si = 2 *  jointLabels[i][i] / ((double) marginalA[i] + marginalB[i] );

        if (siRow == True){
          cout << si << ",";
        } else if (orRow == True){
          cout << overlap << ",";
        } else {
          cout << i << "," << marginalA[i] << "," << marginalB[i] << "," << jointLabels[i][i];
          if (printFalseValues == True){
            cout << "," << marginalA[i] - jointLabels[i][i] << "," << marginalB[i] - jointLabels[i][i];
          }
          cout << "," << overlap << "," << si << endl;
        }

      }
    }

    if (siRow == True || orRow == True){
      cout << endl;
    }
  }

}
