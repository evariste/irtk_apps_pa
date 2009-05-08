#include <irtkImage.h>
#include <irtkTransformation.h>

char *outputFFDName       = NULL;
char *inputFFDName1       = NULL;
char *inputFFDName2       = NULL;

void usage(){
  cerr << "\n  Usage: ffdcompose  [ffdOut] [ffd1] [ffd2]" << endl;
  cerr << "" << endl;
  cerr << "          ffd1                 ffd2" << endl;
  cerr << "    src1 ------> tgt1 = src2  ------> tgt2" << endl;
  cerr << "    " << endl;
  cerr << "          ffdOut" << endl;
  cerr << "    src1 --------> tgt2" << endl;
  cerr << "    " << endl;
  cerr << "    i.e. ffdOut (x) = ffd2 ( ffd1 (x) ) " << endl;
  cerr << "    " << endl;
  cerr << "  ffdOut has the same lattice as ffd2." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkMatrix jac, global1, global2, globalOut;
  int i, j, k, xdim, ydim, zdim, numberOfCPs, count;
  double x, y, z;

  if (argc != 4){
    usage();
  }

  outputFFDName = argv[1];
  argc--;
  argv++;
  inputFFDName1 = argv[1];
  argc--;
  argv++;
  inputFFDName2 = argv[1];
  argc--;
  argv++;


  irtkMultiLevelFreeFormTransformation *mffd_1 = new irtkMultiLevelFreeFormTransformation;
  mffd_1->irtkTransformation::Read(inputFFDName1);

  irtkMultiLevelFreeFormTransformation *mffd_2 = new irtkMultiLevelFreeFormTransformation;
  mffd_2->irtkTransformation::Read(inputFFDName2);

  irtkMultiLevelFreeFormTransformation *mffd_out = new irtkMultiLevelFreeFormTransformation;
  mffd_out->irtkTransformation::Read(inputFFDName2);



  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *affd_out = dynamic_cast<irtkFreeFormTransformation3D *>(mffd_out->GetLocalTransformation(0));
  xdim = affd_out->GetX();
  ydim = affd_out->GetY();
  zdim = affd_out->GetZ();

  numberOfCPs = xdim * ydim * zdim;

  // Space to store the control point displacements.
  double *xdata = new double[numberOfCPs];
  double *ydata = new double[numberOfCPs];
  double *zdata = new double[numberOfCPs];

  // Store the global transformation.
  global1 = mffd_1->irtkAffineTransformation::GetMatrix();
  global2 = mffd_2->irtkAffineTransformation::GetMatrix();
  globalOut = global1 * global2;

  irtkMatrix identity(4, 4);
  identity.Ident();


  mffd_out->irtkAffineTransformation::PutMatrix(identity);

  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){
        affd_out->Put(i, j, k, 0.0, 0.0, 0.0);
      }
    }
  }

  count = 0;

  irtkMatrix cp(4, 1);
  irtkMatrix p(4, 1);

  // Loop over points in source lattice.
  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){

        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target
        // coordinates
        affd_out->LatticeToWorld(x, y, z);

        // Keep a copy of the control point location.
        cp(0, 0) = x;
        cp(1, 0) = y;
        cp(2, 0) = z;
        cp(3, 0) = 1;

        // How is the control point affected by the global
        // transformation for the composed function?
        p = globalOut * cp;

        // Now estimate the full composed transformation for this
        // control point location.

        mffd_2->Transform(x, y, z);
        mffd_1->Transform(x, y, z);

        // Subtract the global transformation to get an estimate
        // for the displacement at the control point.
        xdata[count] = x - p(0, 0);
        ydata[count] = y - p(1, 0);
        zdata[count] = z - p(2, 0);

        ++count;
      }
    }
  }

  // Interpolate the ffd.
  affd_out->Interpolate(xdata, ydata, zdata);

  mffd_out->irtkAffineTransformation::PutMatrix(globalOut);

  mffd_out->irtkTransformation::Write(outputFFDName);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;

}
