#include <irtkImage.h>
#include <irtkTransformation.h>
//#include <nr.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>

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
  bool ok;
  int i, j, k, level;
  bool useActiveOnly = false;
  double x, y, z;
  double sum, max, min, d;
  int count;
  _Status statusX, statusY, statusZ;
  int quiet = false;
  int numberOfValues;
//  float *data;
  gsl_vector *data;

  // Check command line
  if (argc < 2){
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-active") == 0)){
      argc--;
      argv++;
      useActiveOnly = true;
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

//    data = new float[1 + numberOfValues];
//    data[0] = 0;
    data = gsl_vector_alloc(numberOfValues);

    for (i = 0; i < ffd->GetX(); i++){
      for (j = 0; j < ffd->GetY(); j++){
        for (k = 0; k < ffd->GetZ(); k++){

          ffd->GetStatusCP(i, j, k, statusX, statusY, statusZ);

          if (useActiveOnly == true &&
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

//          data[1 + count] = d;
          gsl_vector_set(data, count, d);

          ++count;
        }
      }
    }

//    sort(count, data);
    gsl_sort_vector(data);

    if (count > 0){
      if (quiet){
        cout << level;
        cout << "," << min;
        cout << "," << max;
        cout << "," << sum / ((double) count);
//        cout << "," << data[(1 + count) / 4];
//        cout << "," << data[(1 + count) / 2];
//        cout << "," << data[3 * (1 + count) / 4];
        cout << "," << gsl_vector_get(data, count/4);
        cout << "," << gsl_vector_get(data, count/2);
        cout << "," << gsl_vector_get(data, 3*count/4);
        cout << "," << count << endl;
      } else {
        cout << "Level " << level;
        cout << "\tmin " << min ;
        cout << "\tmax " << max;
        cout << "\tmean: " << sum / ((double) count);
//        cout << "\tQ1 " << data[(1 + count) / 4];
//        cout << "\tmedian " << data[(1 + count) / 2];
//        cout << "\tQ3 " << data[3 * (1 + count) / 4];
        cout << "\tQ1 " << gsl_vector_get(data, count/4);
        cout << "\tmedian " << gsl_vector_get(data, count/2);
        cout << "\tQ3 " << gsl_vector_get(data, 3*count/4);
        cout << "\tn " << count << endl;
      }
    }

//    delete [] data;
    gsl_vector_free(data);

  }

}
