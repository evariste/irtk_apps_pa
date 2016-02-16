#include <irtkTransformation.h>

char *dof_name = NULL;

#define MAX_POINTS 100

void usage()
{
  cerr << "Usage: pointtransform [dof] [x1] [y1] [z1] ... <xn> <yn> <zn>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 5)
    usage();

  int count, i;
  double *x, *y, *z;

  irtkTransformation *transformation = NULL;

  dof_name = argv[1];
  argc--;
  argv++;

  transformation = irtkTransformation::New(dof_name);

  x = new double[MAX_POINTS];
  y = new double[MAX_POINTS];
  z = new double[MAX_POINTS];

  count = 0;

  while (argc > 1 && count < MAX_POINTS){

    x[count] = atof(argv[1]);
    argc--;
    argv++;

    if (argc == 1){
      cerr << "\nPoint components must be a multiple of three.\n" << endl;
      usage();
    }

    y[count] = atof(argv[1]);
    argc--;
    argv++;

    if (argc == 1){
      cerr << "\nPoint components must be a multiple of three.\n" << endl;
      usage();
    }

    z[count] = atof(argv[1]);
    argc--;
    argv++;

    ++count;
  }

  for (i = 0; i < count; ++i){
    cout << x[i] << " " << y[i] << " " << z[i] << " -> ";
    transformation->Transform(x[i], y[i], z[i]);
    cout << x[i] << " " << y[i] << " " << z[i] << endl;
  }

  delete [] x;
  delete [] y;
  delete [] z;

}
