/////////////////////////////////////////////////////////////////////
// Find a topology preservation value for a ffd based on its jacobian
// values at the voxels of a given image.

#include <irtkImage.h>

#include <irtkTransformation.h>

// Default filenames
char *input_name = NULL, *dof_name  = NULL;
char *executableName = NULL;

void usage(char *exeName){
  cerr << "usage : " << exeName << " [image] [ffd] <-options>" << endl;
  cerr << "Find a topology preservation penalty value for a ffd based on" << endl;
  cerr << "its jacobian values at the voxels of a given image." << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-pad value>    Padding value in [image], voxels to ignore." << endl;
  cerr << "<-t value>      Threshold for topology change (default = 3)" << endl;
  cerr << "                Jacobian values <t and >1/t are assumed to be ok and." << endl;  
  cerr << "                contribute zero to the penalty." << endl;
  cerr << "<-q>            Quiet, csv style output only : count,penalty,mean" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  double pad = -1.0 * FLT_MAX;
  double threshold = 3.0;
  bool quiet = false;

  irtkRealImage *input = NULL;

  // Check command line
  executableName = argv[0];
  if (argc < 3){
    usage(executableName);
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  dof_name  = argv[1];
  argc--;
  argv++;


  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dof_name);

  // Read image
  input = new irtkRealImage(input_name);

  while (argc > 1){
    ok = true;
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;
      argv++;
      pad = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-t") == 0)){
      argc--;
      argv++;
      threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage(executableName);
    }
  } 


  int i, j, k;
  irtkPoint p;
  irtkRealPixel *ptr2input;
  double jac;
  double penalty;
  int count;

  penalty = 0.0;
  count   = 0;
  ptr2input = input->GetPointerToVoxels();

  for (k = 0; k < input->GetZ(); k++){
    for (j = 0; j < input->GetY(); j++){
      for (i = 0; i < input->GetX(); i++){

        if (*ptr2input > pad){
          p = irtkPoint(i, j, k);

          // Transform point into world coordinates
          input->ImageToWorld(p);
          jac = mffd->irtkTransformation::Jacobian(p._x, p._y, p._z);

          if (jac < (1.0/threshold) || jac > threshold){
            penalty += log (0.5 * (jac*jac + 1.0 / (jac*jac)));
            ++count;
          }
        }

        ptr2input++;
      }
    }
  }

  if (quiet == true){
    cout << count << "," << penalty << "," << penalty / (double) count << endl;
  } else {
    cout << "count : " << count << ", penalty : " << penalty << ", mean : " << penalty / (double) count << endl;
  }

}
