
#include <irtkImage.h>

#include <irtkTransformation.h>

void usage()
{
  cerr << "Usage: mat2dof [matfile] [doffile] [-invert]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkMatrix matrix;
  irtkAffineTransformation transformation;

  if ((argc != 3) && (argc != 4)){
    usage();
  }

  // Read matrix
  matrix.Read(argv[1]);
  matrix.Print();

  if (argc == 4){
    if (strcmp(argv[3], "-invert") == 0){
      matrix.Invert();
    } else {
      usage();
    }
  }

  // Assign to a transformation. This call should update the parameters.
  transformation.PutMatrix(matrix);

  // Write transformation
  transformation.irtkTransformation::Write(argv[2]);
  transformation.Print();
}
