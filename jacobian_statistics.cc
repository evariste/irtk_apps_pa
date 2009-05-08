///////////////////////////////////////////////////////////////
// jacobian statistics over a labelled image.

#include <irtkImage.h>

#define MAX_LABELS 1000

char *labels_name = NULL, *jac_name = NULL;

void usage(){
  cout << "Usage : jacobian-statistics [labelImage] [jacobianImage] <-options>" << endl;
  cout << "\nFind statistics for jacobians for each label in the segmentation " << endl;
  cout << "[labelImage] . Output is in the format:" << endl;
//   cout << "     label, count, E(jac), E(jac^2) " << endl << endl;;
  cout << "     label, count, E(jac), SD(jac) " << endl << endl;;
  cout << "Options: " << endl;
  cout << "     -log      : Give   E(log(jac)) and SD(log(jac)) instead." << endl;
  cout << "     -noheader : Skip the header." << endl;
  cout << "     -mean     : Only report mean." << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 3)
    usage();

  int pad = 0;

  irtkRealPixel *ptr2jac;
  irtkGreyPixel *ptr2label;

  int i, voxels;
  int labelCounts[MAX_LABELS];
  double jacSums[MAX_LABELS];
  double jacSquaredSums[MAX_LABELS];
  int ok, useLogs = False, useHeader = True, meanOnly = False;
  double logJac;
  double ejac, ejacsq;

  labels_name = argv[1];
  argv++;
  argc--;
  jac_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-log") == 0)){
      argc--;
      argv++;
      useLogs = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-noheader") == 0)){
      argc--;
      argv++;
      useHeader = False;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-mean") == 0)){
      argc--;
      argv++;
      meanOnly = True;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  irtkGreyImage *labelImage = new irtkGreyImage(labels_name);
  irtkRealImage *jacImage = new irtkRealImage(jac_name);

  for (i = 0; i < MAX_LABELS; ++i){
    labelCounts[i] = 0;
    jacSums[i] = jacSquaredSums[i] = 0.0;
  }

  if (labelImage->GetNumberOfVoxels() != jacImage->GetNumberOfVoxels()){
    cerr << "Label image : " << labels_name << " and jacobian image " << jac_name << " must have the same number of voxels." << endl;
    exit(1);
  }

  ptr2label = labelImage->GetPointerToVoxels();
  ptr2jac = jacImage->GetPointerToVoxels();
  voxels = labelImage->GetNumberOfVoxels();

  if(useLogs){

    for (i = 0; i < voxels; ++i){
      if (*ptr2label > pad && *ptr2label < MAX_LABELS && *ptr2jac > 0){
        ++labelCounts[*ptr2label];
        logJac = log(*ptr2jac / 100.0);
        jacSums[*ptr2label] += logJac;
        jacSquaredSums[*ptr2label] += logJac * logJac;
      }
      ++ptr2label;
      ++ptr2jac;
    }

  } else {

    for (i = 0; i < voxels; ++i){
      if (*ptr2label > pad && *ptr2label < MAX_LABELS){
        ++labelCounts[*ptr2label];
        jacSums[*ptr2label] += *ptr2jac;
        jacSquaredSums[*ptr2label] += (*ptr2jac) * (*ptr2jac);
      }
      ++ptr2label;
      ++ptr2jac;
    }

  }

  if (useHeader){
    if(useLogs){
      //      cout << "label,count,sum(log(jac)),sum(log(jac)^2),E(log(jac)),E(log(jac)^2)" << endl;
      cout << "label,count,E(log(jac)),SD(log(jac))" << endl;
    } else {
      //      cout << "label,count,sum(jac),sum(jac^2),E(jac),E(jac^2)" << endl;
      cout << "label,count,E(jac),SD(jac)" << endl;
    }
  }

  for (i = 0; i < MAX_LABELS; ++i){
    if (labelCounts[i] > 0){
      cout << i << ",";
      if (meanOnly){
        cout << jacSums[i] / (double) labelCounts[i] << endl;
      } else {
        cout << labelCounts[i] << ",";
//         cout << jacSums[i] << ",";
//         cout << jacSquaredSums[i] << ",";
        // E (Jac)
        ejac = jacSums[i] / (double) labelCounts[i];
        cout << ejac << ",";
        // E (Jac^2)
        ejacsq = jacSquaredSums[i] / (double) labelCounts[i];
//         cout << ejacsq << ",";
        // SD
        cout << sqrt(ejacsq - ejac * ejac) << endl;
      }
    }
  }


}
