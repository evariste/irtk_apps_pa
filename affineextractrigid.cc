////////////////////////////////////////////////////////////////
#include <irtkTransformation.h>

char *dofin_name = NULL, *dofout_name = NULL;

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
  int ok;
  bool blnQuiet = false;

  if (argc < 3)
    usage();

  dofin_name = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      blnQuiet = true;
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
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

  dofout->irtkTransformation::Write(dofout_name);
}
