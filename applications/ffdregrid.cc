
#include <irtkTransformation.h>

// Input transformation
char *dofin_name  = NULL;

char *template_dof_name = NULL;

// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffdregrid [dofin] [templateDof] [dofout]\n" << endl;
  cerr << "Estimate a dof with the required grid (same as templateDof) that has the same effect" << endl;
  cerr << "as dofin.  Affine part copied over to dofout directly." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k;

  double x, y, z;

  // Check command line
  if (argc != 4){
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  template_dof_name = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transformation
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffdIn = new irtkMultiLevelFreeFormTransformation;
  mffdIn->irtkTransformation::Read(dofin_name);

  irtkMultiLevelFreeFormTransformation *mffdTemplate = new irtkMultiLevelFreeFormTransformation;
  mffdTemplate->irtkTransformation::Read(template_dof_name);
  cout << "done" << endl;

  // Default level is the last level
  if ( mffdTemplate->NumberOfLevels() > 1){
    cerr << "Only implemented for single level dofs." << endl;
    usage();
  }

  irtkMatrix globalMatrix(4, 4);
  globalMatrix = mffdIn->irtkAffineTransformation::GetMatrix();
  mffdTemplate->irtkAffineTransformation::PutMatrix(globalMatrix);

  // Extract first transformation level
  irtkFreeFormTransformation3D *ffd =
	  dynamic_cast<irtkFreeFormTransformation3D *>
		  (mffdTemplate->PopLocalTransformation());

  // Space for storing the control point displacements.
  double* xdata = new double[ffd->NumberOfDOFs() / 3];
  double* ydata = new double[ffd->NumberOfDOFs() / 3];
  double* zdata = new double[ffd->NumberOfDOFs() / 3];

  int count = 0;

  for (k = 0; k < ffd->GetZ(); k++){
     for (j = 0; j < ffd->GetY(); j++){
       for (i = 0; i < ffd->GetX(); i++){

        // Get world coordinates.
        x = i;
        y = j;
        z = k;

        ffd->LatticeToWorld(x, y, z);
        mffdIn->LocalDisplacement(x, y, z);

        xdata[count] = x;
        ydata[count] = y;
        zdata[count] = z;

        ++count;
      }
    }
  }

  ffd->Interpolate(xdata, ydata, zdata);

  mffdTemplate->PushLocalTransformation(ffd);

  // Write transformation
  mffdTemplate->irtkTransformation::Write(dofout_name);
}
