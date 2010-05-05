#include <irtkImage.h>

#include <irtkBiasField.h>

#include <irtkBiasCorrection.h>

#include <irtkGaussianBlurring.h>

char *output;

void usage()
{
  cerr << "Usage: intensitymatch [target] [source] [output] <options>" << endl;
  exit(1);
}

void match(irtkGreyImage &target, irtkGreyImage &source, irtkGreyPixel padding)
{
  int i, n;
  irtkGreyPixel *ptr1, *ptr2;
  double a, b, sx, sy, sxx, syy, sxy;

  n   = 0;
  sx  = 0;
  sy  = 0;
  sxx = 0;
  syy = 0;
  sxy = 0;

  ptr1 = target.GetPointerToVoxels();
  ptr2 = source.GetPointerToVoxels();
  
  for (i = 0; i < target.GetNumberOfVoxels(); i++){
    if ((*ptr1 > padding) && (*ptr2 > padding)){
      sx  += *ptr1;
      sxx += *ptr1 * *ptr1;
      sy  += *ptr2;
      syy += *ptr2 * *ptr2;
      sxy += *ptr1 * *ptr2;
      n++;
    }
    ptr1++;
    ptr2++;
  }
  b = (n*sxy - sx*sy) / (n*sxx - sx*sx);
  a = (sy - b*sx) / n;
  cout << "offset = " << a << ", slope = " << b << endl;

  ptr1 = target.GetPointerToVoxels();
  ptr2 = source.GetPointerToVoxels();
  for (i = 0; i < target.GetNumberOfVoxels(); i++){
    if (*ptr1 > padding){
      *ptr2 = round((*ptr2 - a) / b);      
    } else {
      *ptr2 = *ptr1;
    }
    ptr1++;
    ptr2++;
  }
}

void subtract(irtkGreyImage &target, irtkGreyImage &source, irtkGreyPixel padding)
{
  int i, n;
  irtkGreyPixel *ptr1, *ptr2;

  n = 0;
  double diff_mean = 0;
  double diff_std  = 0;

  ptr1 = target.GetPointerToVoxels();
  ptr2 = source.GetPointerToVoxels();
  for (i = 0; i < target.GetNumberOfVoxels(); i++){
    if (*ptr1 > padding){
      diff_mean += *ptr1 - *ptr2;
      n++;
    } 
    ptr1++;
    ptr2++;
  }
  diff_mean /= n;

  ptr1 = target.GetPointerToVoxels();
  ptr2 = source.GetPointerToVoxels();
  for (i = 0; i < target.GetNumberOfVoxels(); i++){
    if (*ptr1 > padding){
      diff_std += pow((*ptr1 - *ptr2) - diff_mean, 2);
    } 
    ptr1++;
    ptr2++;
  }
  diff_std = sqrt(diff_std / n);
  cout << "Mean difference = " << diff_mean << " (std = " << diff_std << ")" << endl;

  ptr1 = target.GetPointerToVoxels();
  ptr2 = source.GetPointerToVoxels();
  for (i = 0; i < target.GetNumberOfVoxels(); i++){
    if (*ptr1 > padding){
      if (fabs((*ptr1 - *ptr2) - diff_mean) > diff_std){
	*ptr2 = padding;
      } 
    }
    ptr1++;
    ptr2++;
  }
}

int main(int argc, char **argv)
{
  double blur;
  int ok, padding;
  
  if (argc < 4) {
    usage();
    exit(1);
  }
  // Input image
  irtkGreyImage target;
  target.Read(argv[1]);
  argc--;
  argv++;

  // Input image
  irtkGreyImage source;
  source.Read(argv[1]);
  argc--;
  argv++;

  // Output file name
  output = argv[1];
  argc--;
  argv++;

  // Default parameters
  padding    = MIN_GREY;
  blur       = 0;

  // Parblur se remaining parameters
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-blur") == 0)){
      argc--;
      argv++;
      blur = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  if (blur > 0){
    irtkGaussianBlurring<irtkGreyPixel> blurring(blur);
    blurring.SetInput (&source);
    blurring.SetOutput(&source);
    blurring.Run();
  }

  irtkGreyImage image3 = source;

  match(target, source, padding);
  subtract(target, source, padding);

  // Create bias field
  irtkBSplineBiasField *biasfield = new irtkBSplineBiasField(target, 40, 40, 40);

  // Create bias correction filter
  irtkBiasCorrection _biascorrection;
  _biascorrection.SetInput(&source, &target);
  _biascorrection.SetOutput(biasfield);
  _biascorrection.SetPadding(padding);
  _biascorrection.Run();

  // Generate bias corrected image for next iteration
  _biascorrection.Apply(image3);
  
  image3.Write(output);
}


