#include <irtkImage.h>
#include <irtkTransformation.h>


// Input transformation
char *dofin_name = NULL;

// Output transformation
char *affine_out_name = NULL;

void usage(char *exe_name)
{
  cerr << "Usage: " << exe_name << " [mffd_in] [areg_out]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  // Check command line
  if (argc != 3){
    usage(argv[0]);
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  affine_out_name = argv[1];
  argc--;
  argv++;

  // Read transformation no. 1
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);

  irtkMatrix m = mffd->GetMatrix();
  irtkAffineTransformation *affineDof = new irtkAffineTransformation;

  affineDof->PutMatrix(m);

  // Write transformation
  affineDof->irtkTransformation::Write(affine_out_name);

  return 0;
}

