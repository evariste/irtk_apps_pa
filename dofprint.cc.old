#include <irtkImage.h>

#include <irtkTransformation.h>

void usage()
{
  cerr << "Usage: dofprint [dof file]" << endl;
  cerr << "Only prints the rigid / affine components." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation *transformation;

  if (argc != 2){
    usage();
  }

  // Read transformation
  transformation = irtkTransformation::New(argv[1]);

  int i;
  for (i = 0; i < 6; i++){
    cout << transformation->Get(i) << " ";
  }

  if (transformation->NumberOfDOFs() > 6){
    for (i = 6; i < 12; i++){
      cout << transformation->Get(i) << " ";
    }
  }

  cout << endl;

}
