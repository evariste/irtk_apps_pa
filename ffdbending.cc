#include <irtkImage.h>

#include <irtkTransformation.h>
#include <nr.h>


// Default filenames
char *target_name = NULL, *dof_name  = NULL, *mask_name = NULL;


void usage()
{
  cerr << "Usage: ffdbending [target] [ffd]\n" << endl;
  cerr << "-padding [value]  Padding value" << endl;
  cerr << "-mask [image]     Mask image" << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, ok, count;
  double x, y, z, bend, sumBend, padding;
  double jac, sumJac, absLogJac, sumAbsLogJac, topPen, sumTopPen;
  float *jacs, *bends, *absLogJacs, *topPens;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse image
  target_name  = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Initialize padding value
  padding = -1.0 * FLT_MAX;

 
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }    
    if ((ok == False) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  irtkRealImage target;
  target.Read(target_name);
  
  irtkGreyImage mask;
  
  // Set up a mask, 
  if (mask_name == NULL){
    mask.Read(target_name);

    irtkGreyPixel *ptr2mask = mask.GetPointerToVoxels();
    irtkRealPixel *ptr2tgt  = target.GetPointerToVoxels();

    for (i = 0; i < target.GetNumberOfVoxels(); i++){
      if (*ptr2tgt > padding)
        *ptr2mask = 1;
      else
        *ptr2mask = 0;
      
      ++ptr2tgt;
      ++ptr2mask;
    }
  } else {
    mask.Read(mask_name);
    if ((mask.GetX() != target.GetX()) ||
        (mask.GetY() != target.GetY()) ||
        (mask.GetZ() != target.GetZ()) ||
        (mask.GetT() != target.GetT())){
      cerr << "ffdbending: Target and mask dimensions mismatch. Exiting." << endl;
      exit(1);
    }
  }
  
  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd;
  mffd = new irtkMultiLevelFreeFormTransformation;
  
  mffd->irtkTransformation::Read(dof_name);

  count = 0;
  for (k = 0; k < target.GetZ(); k++) {
    for (j = 0; j < target.GetY(); j++) {
      for (i = 0; i < target.GetX(); i++) {
        if (mask.Get(i, j, k) > 0) {
          ++count;
        }
      }
    }
  }

  jacs = new float[1 + count];
  bends = new float[1 + count];
  absLogJacs = new float[1 + count];
  topPens = new float[1 + count];
  
  
  sumBend = 0.0;
  sumJac = 0.0;
  sumAbsLogJac = 0.0;
  sumTopPen = 0.0;
  
  count = 0;
  
  for (k = 0; k < target.GetZ(); k++) {
    for (j = 0; j < target.GetY(); j++) {
      for (i = 0; i < target.GetX(); i++) {

        if (mask.Get(i, j, k) > 0) {
          x = i;
          y = j;
          z = k;

          target.ImageToWorld(x, y, z);
          
          bend = mffd->Bending(x, y, z);
          jac = mffd->irtkTransformation::Jacobian(x, y, z);
          absLogJac = fabs(log(jac));
          topPen = log (0.5 * (jac*jac + 1.0 / (jac*jac)));
          
          sumBend += bend;
          sumJac += jac;
          sumAbsLogJac += absLogJac;
          sumTopPen += topPen;

          bends[1 + count] = bend;
          jacs[1 + count] = jac;
          absLogJacs[1 + count] = absLogJac;
          topPens[1 + count] = topPen;
          
          
          ++count;
        }
      }
    }
  }

  sort(count, bends);
  sort(count, jacs);
  sort(count, absLogJacs);
  sort(count, topPens);

  i = 1 + (int) round(0.5 * (count - 1));

  cout << sumBend / (double(count)) << " ";
  cout << sumJac / (double(count)) << " ";
  cout << sumAbsLogJac / (double(count)) << " ";
  cout << sumTopPen / (double(count)) << " ";

  cout << bends[i] << " ";
  cout << jacs[i] << " ";
  cout << absLogJacs[i] << " ";
  cout << topPens[i] << " ";

  cout << count << endl;

}
