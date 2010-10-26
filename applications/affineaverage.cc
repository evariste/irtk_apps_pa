#include <irtkTransformation.h>

char **dofin_names = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "\t Usage: affineaverage [dofout] [N] [dofin1...N] <-includeRigid>\n" << endl;

  cerr << "\t By default, the rigid parameters of the output transformation" << endl;
  cerr << "\t will be the identity and the scale and skew parameters " << endl;
  cerr << "\t are those of the input transformations averaged." << endl;
  cerr << "\t " << endl;

  cerr << "\t By default, an identity transformation from the target " << endl;
  cerr << "\t to itself is assumed present by default.  This means that" << endl;
  cerr << "\t the averaging will be done over N + 1 transformations." << endl;

  cerr << "\t -includeRigid: By default, only the skew and scale paramaters" << endl;
  cerr << "\t  are actually averaged.  If the output should also reflect " << endl;
  cerr << "\t the average translational and rotational parameters, then" << endl;
  cerr << "\t  set this flag." << endl;
  cerr << "\t " << endl;

  cerr << "\t -noID : Do not assume the N+1_th identity transformation is present." << endl;

  exit(1);
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
  int i, inputCount;
  bool ok;
  irtkMatrix *globalMatrices;
  irtkMatrix globalMatrixAv(4, 4);
  bool ignoreRigid = true;
  int useIdentity = 1;

  if (argc < 4){
    usage();
  }

  // Parse arguments.
  dofout_name = argv[1];
  argc--;
  argv++;
  inputCount = atoi(argv[1]);
  argc--;
  argv++;

  dofin_names = new char *[inputCount];
  for (i = 0; i < inputCount; ++i){
    dofin_names[i] = argv[1];
    argc--;
    argv++;
  }

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-includeRigid") == 0)){
      argc--;
      argv++;
      ignoreRigid = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-noID") == 0)){
      argc--;
      argv++;
      useIdentity = 0;
      cout << "Ignoring extra identity matrix." << endl;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read the transformations.
  irtkAffineTransformation *dofout = new irtkAffineTransformation;

  irtkAffineTransformation **dofin = new irtkAffineTransformation *[inputCount];

  if (useIdentity == 1){
    globalMatrices = new irtkMatrix[inputCount + 1];
  } else {
    globalMatrices = new irtkMatrix[inputCount];
  }

  for (i = 0; i < inputCount; ++i){
    dofin[i] = new irtkAffineTransformation;
    dofin[i]->irtkTransformation::Read(dofin_names[i]);
    globalMatrices[i].Initialize(4, 4);
  }

  // The actual averaging.
  for (i = 0; i < inputCount; ++i){
    if (ignoreRigid){
      resetRigidComponents(dofin[i]);
    }
    globalMatrices[i] = dofin[i]->GetMatrix();
  }

  if (useIdentity == 1){
    // Put an identity matrix at the end of the list of global
    // matrices (for the target to itself).
    globalMatrices[inputCount].Initialize(4, 4);
    globalMatrices[inputCount].Ident();
    cout << "Included identity matrix." << endl;
  }

  if (useIdentity == 1)
    inputCount++;

  globalMatrixAv = FrechetMean(globalMatrices, inputCount, 20);

  dofout->PutMatrix(globalMatrixAv);
  dofout->Print();

  dofout->PutMatrix(globalMatrixAv);
  dofout->irtkTransformation::Write(dofout_name);

}
