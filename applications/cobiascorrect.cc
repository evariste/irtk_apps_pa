//#include <irtkImage.h>
//#include <irtkGaussianBlurring.h>
//#include <irtkBiasField.h>
//#include <irtkBiasCorrection2.h>

#include <irtkIntensityMatcher.h>

char *target_name = NULL;
char *source_name = NULL;
char *tissueMask_name = NULL;

char *tgtout_name = NULL;
char *srcout_name = NULL;

char *biasOutputT = NULL;
char *biasOutputS = NULL;
char *dof_name    = NULL;

void usage()
{
  cerr << "Usage: cobiascorrect [target] [source] [tissueMask] [outputTarget] [outputSource] <options>" << endl;
  cerr << ""  << endl;
  cerr << "Apply bias correction to [target] and [source] simultaneously"  << endl;
  cerr << "and write modified files to [outputTarget] and [outputSource]."  << endl;

  cerr << "The source is resampled to the space of the target (using a dof," << endl;
  cerr << "if specified) before a simultaneous bias correction is done.  " << endl;

  cerr << "All output images and bias fields in target space."  << endl;

  cerr << "-padding [val]  : Padding value in target.  " << endl;
  cerr << "-blur [sigma]   : Preprocess by blurring." << endl;

  cerr << "-target_bias [name] : Write target bias field to [name] " << endl;
  cerr << "-source_bias [name] : write source bias field to [name] " << endl;

  cerr << "-dofin [dof]        : Transformation to map source to target. " << endl;
  cerr << "                      Default is the identity."  << endl;
  cerr << "-csfPercentile [val] : Percentile cut-off to remove csf in  " << endl;
  cerr << "                       homogeneous regions (default 5%)." << endl;
  cerr << "-debug" << endl;

  exit(1);
}


int main(int argc, char **argv)
{
  double blur;
  bool ok;
  int padding;
  double csfPercentile = 5;
  bool debug = false;

  if (argc < 5) {
    usage();
    exit(1);
  }

  irtkTransformation *transformation = NULL;
  irtkGreyImage target;
  irtkGreyImage source;
  irtkGreyImage tissueMask;
  irtkGreyImage trSource;

  irtkImageFunction *interpolator = NULL;
  
  irtkBSplineBiasField *biasfield;

  // Input file names.
  target_name = argv[1];
  argc--;
  argv++;

  source_name = argv[1];
  argc--;
  argv++;

  tissueMask_name = argv[1];
  argc--;
  argv++;

  // Output file names.
  tgtout_name = argv[1];
  argc--;
  argv++;

  srcout_name = argv[1];
  argc--;
  argv++;

  // Input images.
  target.Read(target_name);
  source.Read(source_name);
  tissueMask.Read(tissueMask_name);

  // Default parameters
  padding    = MIN_GREY;
  blur       = 0;

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-blur") == 0)){
      argc--;
      argv++;
      blur = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-target_bias") == 0)){
      argc--;
      argv++;
      biasOutputT = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-source_bias") == 0)){
      argc--;
      argv++;
      biasOutputS = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-csfPercentile") == 0)){
      argc--;
      argv++;
      csfPercentile = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      debug = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  //////////////////////////////////////////////////////////////
  // Transformation and optional blurring stuff.

  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }

  interpolator = new irtkLinearInterpolateImageFunction;

  // Create image transformation
  irtkImageTransformation<irtkGreyPixel> *imagetransformation = 
    new irtkImageTransformation<irtkGreyPixel>;

  trSource = target;

  imagetransformation->SetInput (&source, transformation);
  imagetransformation->SetOutput(&trSource);

  imagetransformation->PutTargetPaddingValue(padding);
  imagetransformation->PutInterpolator(interpolator);

  imagetransformation->Run();

  ////////////////////////////////////////////////////////////////
  // The actual matching:

  irtkIntensityMatcher matcher;

  matcher.SetImage1(&target);
  // Regression is first applied to trSource.
  matcher.SetImage2(&trSource);
 
  matcher.SetTissueMask(&tissueMask);
  matcher.SetPadding(padding);
  matcher.SetCSFPercentile(csfPercentile);
  matcher.SetFieldSpacing(40);
  matcher.SetBlur(blur);
  matcher.SetDebug(debug);
  matcher.Run();


  // Writing:
  target.Write(tgtout_name);
  trSource.Write(srcout_name);
//   source.Write(srcout_name);

  if (biasOutputS != NULL){
    biasfield = matcher.GetBiasField2();
    biasfield->Write(biasOutputS);
  }

  if (biasOutputT != NULL){
    biasfield = matcher.GetBiasField1();
    biasfield->Write(biasOutputT);
  }

}


///////////////////////////////////////////////
// Graveyard:


//   // Match the intensities for the source in native space:
//   double alpha, beta, bias;
//   double x, y, z;
//   irtkGreyPixel *ptr2src;

//   // source \approx _beta * target + _alpha
//   alpha = matcher.GetAlpha();
//   beta  = matcher.GetBeta();

//   biasfield = matcher.GetBiasField2();

//   xdim = source.GetX();
//   ydim = source.GetY();
//   zdim = source.GetZ();
//   ptr2src = source.GetPointerToVoxels();

//   for (k = 0; k < zdim; ++k){
//     for (j = 0; j < ydim; ++j){
//       for (i = 0; i < xdim; ++i){
//         x = i;
//         y = j;
//         z = k;

//         // Source world location.
//         source.ImageToWorld(x, y, z);

//         // Corresponding target location.
//         transformation->Inverse(x, y, z);

//         // Bias as a multiplicative factor.
//         bias = exp( biasfield->Bias(x, y, z) );

//         if (preserve == false){
//           *ptr2src = round (((*ptr2src) - alpha) / beta / bias);
//         } else {
//           *ptr2src = round (alpha + ((*ptr2src) - alpha) / bias );
//         }

//         ++ptr2src;
//       }
//     }
//   }
