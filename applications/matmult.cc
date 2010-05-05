#include <irtkImage.h>

char *matA_name = NULL, *matB_name = NULL, *matout_name = NULL;

void usage(){
  cerr << "Usage: matmult [mat_in1] [mat_in2] [matout]" << endl;
  cerr << ""  << endl;
  exit(1);
}

int main(int argc, char **argv){

  if (argc < 4)
    usage();

  matA_name = argv[1];
  argc--;
  argv++;
  matB_name = argv[1];
  argc--;
  argv++;
  matout_name = argv[1];
  argc--;
  argv++;

  irtkMatrix mA, mB;

  mA.Read(matA_name);
  mB.Read(matB_name);

  cout << matA_name << " :" << endl;
  mA.Print();
  cout << matB_name << " :" << endl;
  mB.Print();

  mA *= mB;

  cout << matA_name << " x " << matB_name << " :" << endl;
  mA.Print();

  mA.Write(matout_name);

}
