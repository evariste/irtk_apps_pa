#include <irtkImage.h>
#include <irtkTransformation.h>

// Input transformation
char *dofin_name = NULL;

// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffdscale [dofin1] [scale factor] [dofout]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, level;
  double x, y, z;
  double scale;

  // Check command line
  if (argc != 4){
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  scale       = atof(argv[1]);
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transformation no. 1
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  cout << "done" << endl;

  for (level = 0; level < mffd->NumberOfLevels(); level++){

    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd =
    	dynamic_cast<irtkFreeFormTransformation3D *>
			(mffd->GetLocalTransformation(level));

    for (i = 0; i < ffd->GetX(); i++){
      for (j = 0; j < ffd->GetY(); j++){
	for (k = 0; k < ffd->GetZ(); k++){

	  // Get control point value on first level
	  ffd->Get(i, j, k, x, y, z);

	  // Add levels and store result in first level
	  ffd->Put(i, j, k, scale*x, scale*y, scale*z);

	}
      }
    }
  }

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}
