#include <irtkImage.h>

char *target_name = NULL;
char *source_name = NULL;
char *output_name = NULL;

void usage()
{
  cerr << "regress [target] [source] <options>" << endl;
  cerr << "Find the least squares regression coefficients for the" << endl;
  cerr << "positive intensities in the source and target." << endl;
  cerr << "Images are assumed to have the same dimensions and are in register. " << endl;
  cerr << "The source is taken to be the predictor (X) variable." << endl;
  cerr << "-output : name for an output file containing the source image" << endl;
  cerr << "          modified by the regression function." << endl;
  cerr << "-pad    : padding value in source and target (default = 0)." << endl;
  exit(1);
}

bool equalSizes(irtkGreyImage &tgt, irtkGreyImage &src){
  return (tgt.GetX() == src.GetX()) && (tgt.GetY() == src.GetY()) && (tgt.GetZ() == src.GetZ());
}

int main(int argc, char **argv)
{
  int ok;
  irtkGreyImage tgt, src, mask;
  irtkGreyPixel *pTgt, *pSrc, *pMask;
  int pad = 0, voxels, i;

  double xbar = 0, ybar = 0, ss_yy = 0, ss_xx = 0, ss_xy = 0, b1, b0;
  double n = 0;

  if (argc < 3){
    usage();
  }

  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-pad") == 0)){
      argc--;      argv++;
      pad = atoi(argv[1]);
      argc--;      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-output") == 0)){
      argc--;      argv++;
      output_name = argv[1];
      argc--;      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  tgt.Read(target_name);
  src.Read(source_name);
  mask.Read(source_name);

  if (!equalSizes(tgt, src)){
    cerr << "Images must have equal sizes.";
    exit(1);
  }

  voxels = tgt.GetNumberOfVoxels();
  pTgt   = tgt.GetPointerToVoxels();
  pSrc   = src.GetPointerToVoxels();
  pMask = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if ((*pTgt > pad) && (*pSrc > pad)){
      ybar += *pTgt;
      xbar += *pSrc;
      ++n;
      *pMask = 1;
    } else {
      *pMask = 0;
    }
    ++pTgt;
    ++pSrc;
    ++pMask;
  }

  xbar /= n;
  ybar /= n;
  cout << "n    = " << n << endl;
  cout << "xbar = " << xbar << endl;
  cout << "ybar = " << ybar << endl;

  pTgt   = tgt.GetPointerToVoxels();
  pSrc   = src.GetPointerToVoxels();
  pMask = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pMask > 0){
      ss_xx += (*pSrc - xbar) * (*pSrc - xbar);
      ss_yy += (*pTgt - ybar) * (*pTgt - ybar);
      ss_xy += (*pSrc - xbar) * (*pTgt - ybar);
    }
    ++pTgt;
    ++pSrc;
    ++pMask;
  }

  b1 = ss_xy/ss_xx;
  b0 = ybar - (b1 * xbar);

  cout << "y-intercept = " << b0 << endl;
  cout << "Gradient    = " << b1 << endl;

  if (output_name != NULL){
    pSrc   = src.GetPointerToVoxels();
    for (i = 0; i < voxels; ++i){
      *pSrc = (int) round(b0 + b1 * (*pSrc));
      ++pSrc;
    }
    src.Write(output_name);
  }

}
