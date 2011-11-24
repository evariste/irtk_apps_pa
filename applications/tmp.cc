

#include <irtkTransformation.h>

// Input transformation
char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: " << endl;
  exit(0);
}

//void f(int &x){
//	x = x * 2;
//}

int main(int argc, char **argv)
{

	int i;

	if (argc < 3){
		usage();
	}

	// dof file with represents a velocity field.
  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name  = argv[1];
  argc--;
  argv++;

  irtkTransformation *transform = irtkTransformation::New(dofin_name);

  irtkMultiLevelFreeFormTransformation *mffd =
  		new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *) transform));

  irtkFreeFormTransformation3D *ffd =
  		dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));

  irtkTransformation *transform2 = irtkTransformation::New(dofin_name);

  irtkMultiLevelFreeFormTransformation *mffd2 =
  		new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *) transform2));

  irtkFreeFormTransformation3D *ffd2 =
  		dynamic_cast<irtkFreeFormTransformation3D *>(mffd2->GetLocalTransformation(0));

  ffd->Print();
  cout << ffd << endl << endl;

//  ffd2->Print();
//  cout << ffd2 << endl << endl;

  ffd->ExponentiateEuler();

  for (i = 0; i < ffd2->NumberOfDOFs(); ++i)
  	ffd2->Put(i, -1 * ffd2->Get(i));


//  ffd->Log();
  mffd->irtkTransformation::Write(dofout_name);

  ffd2->ExponentiateEuler();
  mffd2->irtkTransformation::Write("inv.dof");



//  int a[3] = {23, 3, -2};
//  for (i = 0; i < 3; ++i){
//  	f(a[i]);
//  }
//
//  for (i = 0; i < 3; ++i){
//  	cout << a[i] << endl;
//
//  }

}
