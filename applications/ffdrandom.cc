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
  cerr << "Usage: ffdrandom [dofin] [dofout] <-sigma> <-blur> <-nolimits>" << endl;
  cerr << "Create a random ffd with sigma (default 1) using the same control point " << endl;
  cerr << "lattice as [dofin]. Write result to [dofout]" << endl;
  cerr << "" << endl;
  cerr << "Affine component is identity." << endl;
  cerr << "" << endl;
  cerr << "<sigma> value      Override default sigma (1)." << endl;
  cerr << "" << endl;
  cerr << "<-blur>            Components are blurred first." << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "<-nolimits>        Override the default which sets minimum and maximum component" << endl;
  cerr << "                   values for control points to 0.4 * half a spacing.*"  << endl;
  cerr << "" << endl;
  cerr << "<-no_zero_border>  Do not set the components within 2 lattice points of outer boundary to zero" << endl;
  cerr << "                   This is done by default to allow a smooth transition at the edge." << endl;
  cerr << "" << endl;
  cerr << "*   Choi, Y., Lee, S.: Injectivity conditions of 2D and 3D uniform cubic b-spline functions." << endl;
  cerr << "    Graphical Models 62(6) (2000) 411427" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int xdim, ydim, zdim;
  int i, j, k, level;
  bool ok;
  double x, y, z, xspacing, yspacing, zspacing, minspacing;
  double noise_sigma = 1.0;
  double blur_sigma = 0.0;
  bool threeD = false;
  bool setLimits = true;
  bool zero_border = true;

  double fractionLimit = 0.4;

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
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;      argv++;
      noise_sigma = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-blur") == 0)){
      argc--;      argv++;
      blur_sigma = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nolimits") == 0)){
      setLimits = false;
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-no_zero_border") == 0)){
      zero_border = false;
      argc--;      argv++;
      ok = true;
    }
    if (ok == false){
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

    if (zdim > 1)
    	threeD = true;


    if (setLimits == true){
    	minspacing = xspacing;
      if (minspacing > yspacing)
        minspacing = yspacing;
      if (threeD == true && minspacing > zspacing)
        minspacing = zspacing;
      gaussian_noise->SetMinVal(fractionLimit * minspacing * -1.0);
      gaussian_noise->SetMaxVal(fractionLimit * minspacing);
    }

    noise->SetInput(xcomps);
    noise->SetOutput(xcomps);
    noise->Run();
    noise->SetInput(ycomps);
    noise->SetOutput(ycomps);
    noise->Run();
    noise->SetInput(zcomps);
    noise->SetOutput(zcomps);
    noise->Run();

    if (zero_border){
      cout << "Applying zeros at borders." << endl;

    	if (zdim > 4){
    		for (k = 0; k < 2   ; ++k){
    			for (j = 0; j < ydim; ++j){
    				for (i = 0; i < xdim; ++i){
    					xcomps->Put(i, j, k, 0);
    					xcomps->Put(i, j, zdim - 1 - k, 0);
    					ycomps->Put(i, j, k, 0);
    					ycomps->Put(i, j, zdim - 1 - k, 0);
    					zcomps->Put(i, j, k, 0);
    					zcomps->Put(i, j, zdim - 1 - k, 0);
    				}
    			}
    		}
    	} else {
    		cerr << "Warning : zero border requested but one or more of FFD dimensions is below minimum size needed (4)." << endl;
    	}


    	if (ydim > 4){
    		for (k = 0; k < zdim; ++k){
    			for (j = 0; j < 2   ; ++j){
    				for (i = 0; i < xdim; ++i){
    					xcomps->Put(i, j, k, 0);
    					xcomps->Put(i, ydim - 1 - j, k, 0);
    					ycomps->Put(i, j, k, 0);
    					ycomps->Put(i, ydim - 1 - j, k, 0);
    					zcomps->Put(i, j, k, 0);
    					zcomps->Put(i, ydim - 1 - j, k, 0);
    				}
    			}
    		}
    	} else {
    		cerr << "Warning : zero border requested but one or more of FFD dimensions is below minimum size needed (4)." << endl;
    	}


    	if (xdim > 4){
    		for (k = 0; k < zdim; ++k){
    			for (j = 0; j < ydim; ++j){
    				for (i = 0; i < 2   ; ++i){
    					xcomps->Put(i, j, k, 0);
    					xcomps->Put(xdim - 1 - i, j, k, 0);
    					ycomps->Put(i, j, k, 0);
    					ycomps->Put(xdim - 1 - i, j, k, 0);
    					zcomps->Put(i, j, k, 0);
    					zcomps->Put(xdim - 1 - i, j, k, 0);
    				}
    			}
    		}
    	} else {
    		cerr << "Warning : zero border requested but one or more of FFD dimensions is below minimum size needed (4)." << endl;
    	}
    }

    if (blur_sigma > 0.0){
      cout << "Smoothing components with kernel width of " << gaussianBlurring.GetSigma() << endl;
      gaussianBlurring.SetInput (xcomps);
      gaussianBlurring.SetOutput(xcomps);
      gaussianBlurring.Run();
      gaussianBlurring.SetInput (ycomps);
      gaussianBlurring.SetOutput(ycomps);
      gaussianBlurring.Run();
      gaussianBlurring.SetInput (zcomps);
      gaussianBlurring.SetOutput(zcomps);
      gaussianBlurring.Run();
    }

    for (i = 0; i < xdim; i++){
      for (j = 0; j < ydim; j++){
        for (k = 0; k < zdim; k++){

          x = xcomps->Get(i, j, k);
          y = ycomps->Get(i, j, k);
          z = zcomps->Get(i, j, k);

          if (threeD)
          	ffd->Put(i, j, k, x, y, z);
          else
          	ffd->Put(i, j, k, x, y, 0.0);
        }
      }
    }

    delete xcomps;
    delete ycomps;
    delete zcomps;

  }


  mffd->irtkTransformation::Write(dofout_name);
}
