#include <irtkImage.h>
#include <irtkTransformation.h>
#include <nr.h>

// Input transformation
char *dofin_name = NULL;

// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffdmagnitude [dofin] <-options>" << endl;
  cerr << "Find average magnitude of the control points in a dof." << endl;
  cerr << "Options:" << endl;
  cerr << "-active   Only use active control points." << endl;
  cerr << "-q        Quiet output." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i, j, k, level, useActiveOnly = False;
  double x, y, z;
  double sum, max, min, d;
  int count;
  _Status statusX, statusY, statusZ;
  int quiet = False;
  int numberOfValues;
  float *data;

  // Check command line
  if (argc < 2){
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-active") == 0)){
      argc--;
      argv++;
      useActiveOnly = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = True;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read transformation.
  //  cout << "Reading transformation ... "; cout.flush();
  irtkTransformation *transform = irtkTransformation::New(dofin_name);
  irtkMultiLevelFreeFormTransformation *mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);
  //  cout << "done" << endl;

  max = -1.0 * FLT_MAX;
  min = FLT_MAX;

  for (level = 0; level < mffd->NumberOfLevels(); level++){

    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd =
    	dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(level));
    sum = 0.0;
    count = 0;

    numberOfValues = ffd->NumberOfDOFs() / 3;
    data = new float[1 + numberOfValues];
    data[0] = 0;

    for (i = 0; i < ffd->GetX(); i++){
      for (j = 0; j < ffd->GetY(); j++){
	for (k = 0; k < ffd->GetZ(); k++){

          ffd->GetStatus(i, j, k, statusX, statusY, statusZ);

          if (useActiveOnly == True &&
              statusX == _Passive && statusY == _Passive && statusZ == _Passive){
            continue;
          }

	  // Get control point value on first level
	  ffd->Get(i, j, k, x, y, z);
          d = sqrt(x*x + y*y + z*z);

          if (max < d)
            max = d;
          if (min > d)
            min = d;

          sum += d;

          data[1 + count] = d;

          ++count;
	}
      }
    }

    sort(count, data);

    if (count > 0){
      if (quiet){
        cout << level;
        cout << "," << min;
        cout << "," << max;
        cout << "," << sum / ((double) count);
        cout << "," << data[(1 + count) / 4];
        cout << "," << data[(1 + count) / 2];
        cout << "," << data[3 * (1 + count) / 4];
        cout << "," << count << endl;
      } else {
        cout << "Level " << level;
        cout << "\tmin " << min ;
        cout << "\tmax " << max;
        cout << "\tmean: " << sum / ((double) count);
        cout << "\tQ1 " << data[(1 + count) / 4];
        cout << "\tmedian " << data[(1 + count) / 2];
        cout << "\tQ3 " << data[3 * (1 + count) / 4];
        cout << "\tn " << count << endl;
      }
    }

    delete [] data;

  }

}
