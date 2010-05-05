#include <irtkImage.h>
#include <irtkNoise.h>
#include <irtkGaussianBlurring.h>
#include <irtkTransformation.h>

// Input transformation
char *dofin_name = NULL;
// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffdrandom [dofin] [dofout] <-sigma> <-blur> <-3d> <-nolimits>" << endl;
  cerr << "Create a random ffd with sigma (default 1) at the same control points." << endl;
  cerr << "Affine component is identity.  If blur specified then components are blurred first." << endl;
  cerr << "Default is for a 2D image unless -3d specified, then z components also included." << endl;
  cerr << "Default is to set min and max values for random cp components as +/- half a spacing."  << endl;
  cerr << "This is overridden if -nolimits is specified."  << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int xdim, ydim, zdim;
  int i, j, k, level, ok;
  double x, y, z, xspacing, yspacing, zspacing, minspacing;
  double noise_sigma = 1.0;
  double blur_sigma = 0.0;
  int threeD = False;
  int setLimits = True;

  irtkRealImage *xcomps, *ycomps, *zcomps;

  if (argc < 3)
    usage();

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;      argv++;
      noise_sigma = atof(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-blur") == 0)){
      argc--;      argv++;
      blur_sigma = atof(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-3d") == 0)){
      threeD = True;
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-nolimits") == 0)){
      setLimits = False;
      argc--;      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkTransformation *transform = irtkTransformation::New(dofin_name);

  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
  delete transform;

  mffd->PutTranslationX(0.0);
  mffd->PutTranslationY(0.0);
  mffd->PutTranslationZ(0.0);
  mffd->PutRotationX(0.0);
  mffd->PutRotationY(0.0);
  mffd->PutRotationZ(0.0);
  mffd->PutScaleX(100.0);
  mffd->PutScaleY(100.0);
  mffd->PutScaleZ(100.0);
  mffd->PutShearXY(0.0);
  mffd->PutShearYZ(0.0);
  mffd->PutShearXZ(0.0);

  irtkNoise<irtkRealPixel> *noise = NULL;
  irtkGaussianNoise<irtkRealPixel> *gaussian_noise = NULL;
  noise = gaussian_noise = new irtkGaussianNoise<irtkRealPixel>;
  gaussian_noise->SetMean(0.0);
  gaussian_noise->SetSigma(noise_sigma);

  irtkGaussianBlurring<irtkRealPixel> gaussianBlurring(blur_sigma);

  irtkImageAttributes attr;

  for (level = 0; level < mffd->NumberOfLevels(); level++){
    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd =
    	dynamic_cast<irtkFreeFormTransformation3D *> (mffd->GetLocalTransformation(level));
    ffd->GetSpacing(xspacing, yspacing, zspacing);
    xdim = ffd->GetX();
    ydim = ffd->GetY();
    zdim = ffd->GetZ();

    attr._x = xdim;
    attr._y = ydim;
    attr._z = zdim;
    attr._dx = xspacing;
    attr._dy = yspacing;
    attr._dz = zspacing;

    xcomps = new irtkRealImage(attr);
    ycomps = new irtkRealImage(attr);
    zcomps = new irtkRealImage(attr);

//    xcomps = new irtkRealImage(xdim, ydim, zdim, xspacing, yspacing, zspacing);
//    ycomps = new irtkRealImage(xdim, ydim, zdim, xspacing, yspacing, zspacing);
//    zcomps = new irtkRealImage(xdim, ydim, zdim, xspacing, yspacing, zspacing);

    if (setLimits == True){
      minspacing = xspacing;
      if (minspacing > yspacing)
        minspacing = yspacing;
      if (threeD == True && minspacing > zspacing)
        minspacing = zspacing;
      gaussian_noise->SetMinVal(-0.5*minspacing);
      gaussian_noise->SetMaxVal(0.5*minspacing);
    }

    noise->SetInput(xcomps);
    noise->SetOutput(xcomps);
    noise->Run();
    noise->SetInput(ycomps);
    noise->SetOutput(ycomps);
    noise->Run();
    if (threeD == True){
      noise->SetInput(zcomps);
      noise->SetOutput(zcomps);
      noise->Run();
    }

    if (blur_sigma > 0.0){
      gaussianBlurring.SetInput (xcomps);
      gaussianBlurring.SetOutput(xcomps);
      gaussianBlurring.Run();
      gaussianBlurring.SetInput (ycomps);
      gaussianBlurring.SetOutput(ycomps);
      gaussianBlurring.Run();
      if (threeD == True){
        gaussianBlurring.SetInput (zcomps);
        gaussianBlurring.SetOutput(zcomps);
        gaussianBlurring.Run();
      }
    }
//     xcomps->Write("xcomps.gipl");
//     ycomps->Write("ycomps.gipl");
//     zcomps->Write("zcomps.gipl");

    for (i = 0; i < xdim; i++){
      for (j = 0; j < ydim; j++){
        for (k = 0; k < zdim; k++){

          x = xcomps->Get(i, j, k);
          y = ycomps->Get(i, j, k);
          z = zcomps->Get(i, j, k);

          ffd->Put(i, j, k, x, y, z);
        }
      }
    }

    delete xcomps;
    delete ycomps;
    delete zcomps;

  }


  mffd->irtkTransformation::Write(dofout_name);
}
