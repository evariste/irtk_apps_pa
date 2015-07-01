#include <irtkTransformation.h>

char *dofin_name1 = NULL;
char *dofin_name2 = NULL;
char *affine_out_name = NULL;
      
void usage()
{
  cerr << "Usage: affinecompose [doffile1] [doffile2] [doffileOut]" << endl;
  cerr << "     doffile1              doffile2" << endl;
  cerr << "src1 <------- tgt1 = src2  <------- tgt2" << endl;
  cerr << "" << endl;
  cerr << "     doffileOut" << endl;
  cerr << "src1 <--------- tgt2" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  irtkMatrix matrix1, matrix2, matrixOut;
  irtkAffineTransformation transformation1, transformation2;

  if (argc < 2){
    usage();
  }

  // Read transformations
  dofin_name1 = argv[1];
  argc--;
  argv++;
  dofin_name2 = argv[1];
  argc--;
  argv++;
  affine_out_name = argv[1];
  argc--;
  argv++;

  transformation1.irtkTransformation::Read(dofin_name1);
  transformation2.irtkTransformation::Read(dofin_name2);

  // Convert to matrices
  matrix1 = transformation1.GetMatrix();
  matrix2 = transformation2.GetMatrix();

  matrixOut = matrix1 * matrix2;
  transformation1.PutMatrix(matrixOut);

  transformation1.irtkTransformation::Write(affine_out_name);

}
