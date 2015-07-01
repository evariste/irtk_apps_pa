////////////////////////////////////////////////////////////////
#include <irtkTransformation.h>

char *dofin_name = NULL, *affine_out_name = NULL;

void usage()
{
  cerr << "Usage: affineextractrigid [dof_in] [dof_out] <-q>" << endl;
  cerr << "-q : quiet." << endl;
  cerr << ""  << endl;
  exit(1);
}

void myPrint(irtkAffineTransformation *t)
{
  cout << "Number of DOFs: " << t->NumberOfDOFs() << endl;
  cout << "\t\t -tx " << t->GetTranslationX();
  cout << " -ty " << t->GetTranslationY();
  cout << " -tz " << t->GetTranslationZ();
  cout << " -rx " << t->GetRotationX();
  cout << " -ry " << t->GetRotationY();
  cout << " -rz " << t->GetRotationZ();
  cout << endl;

  if (t->NumberOfDOFs() > 6){
    cout << "\t\t -sx " << t->GetScaleX();
    cout << " -sy " << t->GetScaleY();
    cout << " -sz " << t->GetScaleZ();
    cout << " -sxy " << t->GetShearXY();
    cout << " -sxz " << t->GetShearXZ();
    cout << " -syz " << t->GetShearYZ();
    cout << endl;
  }
}

void myPrint(irtkRigidTransformation *t)
{
  cout << "Number of DOFs: " << t->NumberOfDOFs() << endl;
  cout << "\t\t -tx " << t->GetTranslationX();
  cout << " -ty " << t->GetTranslationY();
  cout << " -tz " << t->GetTranslationZ();
  cout << " -rx " << t->GetRotationX();
  cout << " -ry " << t->GetRotationY();
  cout << " -rz " << t->GetRotationZ();
  cout << endl;
}

int main(int argc, char **argv)
{
  bool ok;
  bool blnQuiet = false;

  if (argc < 3)
    usage();

  dofin_name = argv[1];
  argc--;
  argv++;
  affine_out_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      blnQuiet = true;
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkTransformation *transform;

  irtkAffineTransformation *dofin;

  transform = irtkTransformation::New(dofin_name);

  dofin  = new irtkAffineTransformation(*((irtkAffineTransformation *)transform));

  irtkRigidTransformation *dofout;

  dofout = new irtkRigidTransformation;

  for (int i = 0; i < 6; i++){
    dofout->Put(i, dofin->Get(i));
  }

  if (!blnQuiet){
    cout << "Dofin  : " << endl;
    myPrint(dofin);
    cout << "Dofout  : " << endl;
    myPrint(dofout);
  }

  dofout->irtkTransformation::Write(affine_out_name);
}
