#include <irtkTransformation.h>

char **dofin_names = NULL, *dofout = NULL;

void usage()
{
  cerr << "\t Usage: ffdaverage [dofout] [N] [dofin1...N] <-noaffine> \n" << endl;

  cerr << "\t The input mffds are assumed to relate to the same target image" << endl;
  cerr << "\t and have the same lattice dimensions and a single level. " << endl;
  cerr << "\t The output mffd will also point to this target and it will " << endl;
  cerr << "\t have a single level with the same control point structure as the " << endl;
  cerr << "\t input transformations." << endl;
  cerr << "\t The rigid parameters of the output mffd will be the identity and " << endl;
  cerr << "\t the scale and skew parameters are those of the input mffds averaged." << endl;
  cerr << "\t " << endl;
  cerr << "\t The input transformation control points are averaged linearly " << endl;
  cerr << "\t and assigned to each output control point.  Affine correction is first " << endl;
  cerr << "\t done on the input transformations' local control point values." << endl ;
  cerr << "\t An identity mffd from the target to itself is assumed present." << endl;
  cerr << "\t -noaffine     = Don't combine with affine average, i.e. resulting " << endl;
  cerr << "\t                 transformation will have identity global component." << endl;
  exit(1);
}

void checkLevelsAndDimensions(irtkMultiLevelFreeFormTransformation **mffds, int inputCount){
  int i, x, y, z;

  x = mffds[0]->GetLocalTransformation(0)->GetX();
  y = mffds[0]->GetLocalTransformation(0)->GetY();
  z = mffds[0]->GetLocalTransformation(0)->GetZ();

  for (i = 0; i < inputCount; ++i){
    if (mffds[i]->NumberOfLevels() > 1){
      cerr << "All input transformations must have the same number of levels." << endl;
      exit(1);
    }
    if (x != mffds[i]->GetLocalTransformation(0)->GetX() ||
        y != mffds[i]->GetLocalTransformation(0)->GetY() ||
        z != mffds[i]->GetLocalTransformation(0)->GetZ()){
      cerr << "All input transformations must have the same lattice dimensions." << endl;
      exit(1);
    }
  }
}

void resetRigidComponents(irtkAffineTransformation *a){
  a->PutTranslationX(0);
  a->PutTranslationY(0);
  a->PutTranslationZ(0);
  a->PutRotationX(0);
  a->PutRotationY(0);
  a->PutRotationZ(0);
}

int main(int argc, char **argv)
{
  bool ok;
  int i, j, k, n, xdim, ydim, zdim, inputCount;
  double x, y, z, xAffCorr, yAffCorr, zAffCorr;
  double dx, dy, dz;
  double xAv, yAv, zAv, xFinal, yFinal, zFinal;
  int combineaffine = true;

  irtkMatrix *globalMatrices;
  irtkMatrix *localMatrices;

  irtkMatrix globalMatrixAv(4, 4);
  irtkMatrix sumLogs(4, 4);
  irtkMatrix jac(3, 3);
  irtkMatrix globalAvJac(3, 3);

  if (argc < 4){
    usage();
  }

  // Parse arguments.
  dofout = argv[1];
  argc--;
  argv++;
  inputCount = atoi(argv[1]);
  argc--;
  argv++;

  if (inputCount < 1){
    cerr << "checkLevelsAndDimensions: Must have at least one input transformation." << endl;
    exit(1);
  }


  // Read input transformation names.
  dofin_names = new char *[inputCount];
  for (i = 0; i < inputCount; ++i){
    dofin_names[i] = argv[1];
    argc--;
    argv++;
  }

  // Parse options.
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-noaffine") == 0)){
      argc--;
      argv++;
      combineaffine = false;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read the mffds.
  // The output mffd is a copy of the first input (for now).
  irtkTransformation *transform;
  irtkMultiLevelFreeFormTransformation *mffdOut;
  irtkMultiLevelFreeFormTransformation **mffds;

  transform = irtkTransformation::New(dofin_names[0]);
  mffdOut   = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);

  mffds          = new irtkMultiLevelFreeFormTransformation *[inputCount];
  globalMatrices = new irtkMatrix[inputCount + 1];

  for (i = 0; i < inputCount; ++i){
    transform = irtkTransformation::New(dofin_names[i]);
    mffds[i]  = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);
    globalMatrices[i].Initialize(4, 4);
  }

  checkLevelsAndDimensions(mffds, inputCount);

  ////////////////////////////////////////////////////
  // Affine part:

  globalMatrixAv.Ident();

  // Is affine portion to be averaged?
  if (combineaffine == true){

    // Average the affine parts separately from the displacement fields.
    // Only the scales and skew are actually averaged.
    for (i = 0; i < inputCount; ++i){

      resetRigidComponents(mffds[i]);
      globalMatrices[i] = mffds[i]->irtkAffineTransformation::GetMatrix();
    }

    // Put the identity matrix at the end.
    globalMatrices[inputCount].Initialize(4, 4);
    globalMatrices[inputCount].Ident();

    globalMatrixAv = FrechetMean(globalMatrices, 1 + inputCount, 20);

  }

  mffdOut->irtkAffineTransformation::PutMatrix(globalMatrixAv);

  ////////////////////////////////////////////////////
  // Now tackle the displacement fields.

  // Re-read the original global matrices.
  for (i = 0; i < inputCount; ++i){
    transform = irtkTransformation::New(dofin_names[i]);
    mffds[i]  = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);
    globalMatrices[i] = mffds[i]->irtkAffineTransformation::GetMatrix();
  }


  // Space for storing the output control points.
  irtkFreeFormTransformation3D *ffdOut =
	  dynamic_cast<irtkFreeFormTransformation3D *>(mffdOut->GetLocalTransformation(0));

  // Space for storing the control point displacements.
  double* xdata = new double[ffdOut->NumberOfDOFs()/3];
  double* ydata = new double[ffdOut->NumberOfDOFs()/3];
  double* zdata = new double[ffdOut->NumberOfDOFs()/3];

  xdim = ffdOut->GetX();
  ydim = ffdOut->GetY();
  zdim = ffdOut->GetZ();


  // Correct the control point values in each of the input FFDs with
  // respect to their corresponding affine components.
  for (n = 0; n < inputCount; ++n){

    // Correcting the cp values requires the affine part to be S->T.
    globalMatrices[n].Invert();

    for (k = 0; k < zdim; ++k){
      for (j = 0; j < ydim; ++j){
        for (i = 0; i < xdim; ++i){
          mffds[n]->GetLocalTransformation(0)->Get(i, j, k, x, y, z);

          xAffCorr = globalMatrices[n](0, 0) * x + globalMatrices[n](0, 1) * y + globalMatrices[n](0, 2) * z;
          yAffCorr = globalMatrices[n](1, 0) * x + globalMatrices[n](1, 1) * y + globalMatrices[n](1, 2) * z;
          zAffCorr = globalMatrices[n](2, 0) * x + globalMatrices[n](2, 1) * y + globalMatrices[n](2, 2) * z;

          mffds[n]->GetLocalTransformation(0)->Put(i, j, k, xAffCorr, yAffCorr, zAffCorr);
        }
      }
    }

  }

  // Copy the average jacobian section for correcting the CP values in the
  // output transformation.
  for (i = 0; i < 3; ++i){
    for (j = 0; j < 3; ++j){
      globalAvJac(i, j) = globalMatrixAv(i, j);
    }
  }

  irtkMatrix temp(4, 4);
  irtkMatrix affEstimate(4, 4);
  double wx, wy, wz;
  bool folding = false;
  int foldingCount = 0;
  int count = 0;

  ffdOut->GetSpacing(dx, dy, dz);

  localMatrices = new irtkMatrix[1 + inputCount];

  // Initialize the local matrices.
  for(i = 0; i < inputCount; i++){
    localMatrices[i].Initialize(4, 4);
  }
  // Put the identity matrix at the end.
  localMatrices[inputCount].Initialize(4, 4);
  localMatrices[inputCount].Ident();

  // Now can average the displacement fields.
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){

        // Find current world coordinates.
        wx = i;    wy = j;    wz = k;
        ffdOut->LatticeToWorld(wx, wy, wz);

        for (n = 0; n < inputCount; ++n){

          affEstimate.Ident();
          x = wx;    y = wy;   z = wz;

          mffds[n]->GetLocalTransformation(0)->LocalJacobian(x, y, z, jac);

          if (jac.Det() < 0.001){
            jac.Ident();
            folding = true;
          }

          // jac is 3x3 copy it into the top left of affEstimate.
          affEstimate(jac, 0, 0);

          // Put the translation components into affEstimate.
          mffds[n]->GetLocalTransformation(0)->LocalDisplacement(x, y, z);
          affEstimate(0, 3) = x;
          affEstimate(1, 3) = y;
          affEstimate(2, 3) = z;

          localMatrices[n] = affEstimate;
        }

        temp = FrechetMean(localMatrices, 1 + inputCount, 20);

        xAv = temp(0, 3);
        yAv = temp(1, 3);
        zAv = temp(2, 3);

        if (combineaffine == true){
          // Now introduce the affine component.
          xFinal = globalAvJac(0, 0) * xAv + globalAvJac(0, 1) * yAv + globalAvJac(0, 2) * zAv;
          yFinal = globalAvJac(1, 0) * xAv + globalAvJac(1, 1) * yAv + globalAvJac(1, 2) * zAv;
          zFinal = globalAvJac(2, 0) * xAv + globalAvJac(2, 1) * yAv + globalAvJac(2, 2) * zAv;
        } else {
          xFinal = xAv;
          yFinal = yAv;
          zFinal = zAv;
        }

        ffdOut->Put(i, j, k, xFinal, yFinal, zFinal);

        if (folding){
          ++foldingCount;
          // Reset.
          folding = false;
        }
        ++count;
      }
    }
  }

  // Some reporting.
  if (combineaffine == true){
    cout << "Affine component combined with ffd." << endl;
  } else {
    cout << "ffd has identity affine component." << endl;
  }
  cout << "Control points: " << count << endl;
  cout << "Singularities:  " << foldingCount << endl;

  mffdOut->irtkTransformation::Write(dofout);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;
}
