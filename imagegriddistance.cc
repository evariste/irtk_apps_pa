#include <irtkImage.h>
#include <irtkTransformation.h>

char *target_name = NULL, *source_name = NULL;
char *dof_name = NULL;

void usage()
{
  cerr << "\timagegriddistance [target] [source] <-dofin dof>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  int i, j, k, ok;
  double x, y, z;
  double l, m, n, s, t, u;
  double sum, sumx, sumy, sumz;

  irtkGreyImage *target, *source;

  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  target = new irtkGreyImage(target_name);
  source = new irtkGreyImage(source_name);

  irtkTransformation *transformation = NULL;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  int xdim = target->GetX();
  int ydim = target->GetY();
  int zdim = target->GetZ();


  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }

  sum = sumx = sumy = sumz = 0.0;

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){

        x = i;
        y = j;
        z = k;
        target->ImageToWorld(x, y, z);
        transformation->Transform(x, y, z);
        source->WorldToImage(x, y, z);

        l = (int)round(x);
        m = (int)round(y);
        n = (int)round(z);

        s = x - l;
        t = y - m;
        u = z - n;

        sum += s*s + t*t + u*u;
        sumx += s*s;
        sumy += t*t;
        sumz += u*u;

      }
    }
  }

  sum /= (double) target->GetNumberOfVoxels();
  sumx /= (double) target->GetNumberOfVoxels();
  sumy /= (double) target->GetNumberOfVoxels();
  sumz /= (double) target->GetNumberOfVoxels();

  sum = sqrt(sum);
  sumx = sqrt(sumx);
  sumy = sqrt(sumy);
  sumz = sqrt(sumz);

  cout << sumx << " " << sumy << " "  << sumz << " " << sum << endl;

}



