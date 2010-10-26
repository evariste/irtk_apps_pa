////////////////////////////////////////////////////////////////
#include <irtkTransformation.h>

char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: affineresetrigid [dof_in] [dof_out] <-q>" << endl;
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

int main(int argc, char **argv)
{
  bool ok;
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


  // Just in case the scales etc have zero or garbage in them ...
  if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0){
    if (!blnQuiet){
      cout << "Converting from rigid to affine." << endl;
    }
    dofin->Put(6, 100);
    dofin->Put(7, 100);
    dofin->Put(8, 100);

    dofin->Put(9, 0);
    dofin->Put(10, 0);
    dofin->Put(11, 0);
  }

  if (!blnQuiet){
    cout << "Dofin  : " << endl;
    myPrint(dofin);
  }

  dofin->PutTranslationX(0);
  dofin->PutTranslationY(0);
  dofin->PutTranslationZ(0);
  dofin->PutRotationX(0);
  dofin->PutRotationY(0);
  dofin->PutRotationZ(0);

  if (!blnQuiet){
    cout << "Dofout : " << endl;
    myPrint(dofin);
  }

  dofin->irtkTransformation::Write(dofout_name);
}
