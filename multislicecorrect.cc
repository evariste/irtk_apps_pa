#include <irtkImage.h>
#include <irtkResampling.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

char *_template_name  = NULL;
char *_output_name  = NULL;
char **_input_name  = NULL;
char **_dof_name    = NULL;
char *_evalTarget_name   = NULL;

irtkRealImage **_dx   = NULL;
irtkRealImage **_coeffs  = NULL;
irtkRealImage **_weights = NULL;

irtkRealImage *_temp = NULL;
irtkGreyImage *_temp2 = NULL;

irtkRealImage **_input      = NULL;
irtkRealImage *_template    = NULL;
irtkRealImage *_evalTarget  = NULL;

irtkTransformation **_transf = NULL;

int _inputCount = 0;

#define LOOKUPTABLESIZE 1000

#define DBL_LUTSIZE (double)(LOOKUPTABLESIZE-1)

double LookupTable  [LOOKUPTABLESIZE][4];

#define MAX_IMAGES 150

inline double Bsp(int i, double t)
{
  switch (i) {
  case 0:
    return (1-t)*(1-t)*(1-t)/6.0;
  case 1:
    return (3*t*t*t - 6*t*t + 4)/6.0;
  case 2:
    return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
  case 3:
    return (t*t*t)/6.0;
  }
  return 0;
}

void initialiseLookupTable()
{
  int i;

  for (i = 0; i < LOOKUPTABLESIZE; ++i){
    LookupTable[i][0] = Bsp(0, i/DBL_LUTSIZE);
    LookupTable[i][1] = Bsp(1, i/DBL_LUTSIZE);
    LookupTable[i][2] = Bsp(2, i/DBL_LUTSIZE);
    LookupTable[i][3] = Bsp(3, i/DBL_LUTSIZE);
  }
}

void clearRealImage(irtkRealImage *img)
{
  int n = img->GetNumberOfVoxels();
  irtkRealPixel *ptr2pix = img->GetPointerToVoxels();
  memset((void*) ptr2pix, 0, n * sizeof(irtkRealPixel));
}

void clearGreyImage(irtkGreyImage *img)
{
  int n = img->GetNumberOfVoxels();
  irtkGreyPixel *ptr2pix = img->GetPointerToVoxels();
  memset((void*) ptr2pix, 0, n * sizeof(irtkGreyPixel));
}

irtkRealImage * zeroPad(irtkRealImage *in, int count)
{
  int i, j, k;
  double xspacing, yspacing, zspacing;
  double xaxis[3], yaxis[3], zaxis[3];
  irtkPoint origin;
  int xinc, yinc, zinc;
  int xdim, ydim, zdim;
  in->GetPixelSize(&xspacing, &yspacing, &zspacing);
  origin = in->GetOrigin();
  in->GetOrientation(xaxis, yaxis, zaxis);

  xdim = in->GetX();
  ydim = in->GetY();
  zdim = in->GetZ();

  xinc = yinc = zinc = 2 * count;
  // Create new image with required dimensions.
  irtkRealImage *out;

  irtkImageAttributes attr;
  attr._x = xdim + xinc;
  attr._y = ydim + yinc;
  attr._z = zdim + zinc;

  attr._dx = xspacing;
  attr._dy = yspacing;
  attr._dz = zspacing;

  attr._xorigin = origin._x;
  attr._yorigin = origin._y;
  attr._zorigin = origin._z;

  for (i = 0; i < 3; ++i){
    attr._xaxis[i] = xaxis[i];
    attr._yaxis[i] = yaxis[i];
    attr._zaxis[i] = zaxis[i];
  }
  out = new irtkRealImage(attr);
//  out = new irtkRealImage(xdim + xinc, ydim + yinc, zdim + zinc,
//                         xspacing, yspacing, zspacing,
//                         origin, xaxis, yaxis, zaxis);

  // Copy data:
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        out->Put(count + i, count + j, count + k, in->Get(i, j, k));
      }
    }
  }

  return out;
}


// Append, rows, cols, etc. so that the resulting image has
// dimensions that are all 1 mod (16).  This allows the image to
// undergo up to 'levels' downsamplings and subsequent FFD style
// subdivisions.
irtkRealImage * growImage(irtkRealImage *in, int levels)
{
  int i, j, k;
  double xspacing, yspacing, zspacing;
  double xaxis[3], yaxis[3], zaxis[3];
  irtkPoint origin;
  int xinc, yinc, zinc;
  int xdim, ydim, zdim;
  in->GetPixelSize(&xspacing, &yspacing, &zspacing);
  origin = in->GetOrigin();
  in->GetOrientation(xaxis, yaxis, zaxis);

  int requiredBase = (int) round(pow(2.0, levels));

  xdim = in->GetX();
  ydim = in->GetY();
  zdim = in->GetZ();

  // How much must each dimension be incremented by?
  xinc = (requiredBase + 1 - (xdim % requiredBase)) % requiredBase;
  yinc = (requiredBase + 1 - (ydim % requiredBase)) % requiredBase;
  zinc = (requiredBase + 1 - (zdim % requiredBase)) % requiredBase;

  // Get the current world to image matrix.
  irtkMatrix w2i = in->GetWorldToImageMatrix();

  irtkMatrix tmp1(4, 4);
  irtkMatrix tmp2(4, 4);

  // Matrix that gives the voxel lattice an origin at the corner.
  tmp1.Ident();
  tmp1(0, 3) = (xdim - 1) / 2.0;
  tmp1(1, 3) = (ydim - 1) / 2.0;
  tmp1(2, 3) = (zdim - 1) / 2.0;

  // Matrix with the displacement between the centre of the
  // lattice and the origin.
  tmp2.Ident();
  tmp2(0, 3) = - origin._x;
  tmp2(1, 3) = - origin._y;
  tmp2(2, 3) = - origin._z;

  // Calculate the scales and rotations part of w2i matrix.
  irtkMatrix SR(4, 4);

  tmp1.Invert();
  tmp2.Invert();
  SR = tmp1 * (w2i * tmp2);

  // Voxel lattice origin shift matrix for new dimensions.
  tmp1.Ident();
  tmp1(0, 3) = (xdim + xinc - 1) / 2.0;
  tmp1(1, 3) = (ydim + yinc - 1) / 2.0;
  tmp1(2, 3) = (zdim + zinc - 1) / 2.0;

  // Calculate world origin displacement for new size lattice.
  tmp1.Invert();
  SR.Invert();
  tmp2 = SR * tmp1 * w2i;

  origin._x = - tmp2(0, 3);
  origin._y = - tmp2(1, 3);
  origin._z = - tmp2(2, 3);

  // Create new image with required dimensions.
  irtkRealImage *out;
  irtkImageAttributes attr;
  attr._x = xdim + xinc;
  attr._y = ydim + yinc;
  attr._z = zdim + zinc;

  attr._dx = xspacing;
  attr._dy = yspacing;
  attr._dz = zspacing;

  attr._xorigin = origin._x;
  attr._yorigin = origin._y;
  attr._zorigin = origin._z;

  for (i = 0; i < 3; ++i){
    attr._xaxis[i] = xaxis[i];
    attr._yaxis[i] = yaxis[i];
    attr._zaxis[i] = zaxis[i];
  }
  out = new irtkRealImage(attr);
//  out = new irtkRealImage(xdim + xinc, ydim + yinc, zdim + zinc,
//                         xspacing, yspacing, zspacing,
//                         origin, xaxis, yaxis, zaxis);

  // Testing stuff:
  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        out->Put(i, j, k, in->Get(i, j, k));
      }
    }
  }

  return out;
}

// Downsample in way such that FFD style subdivision applied to
// the output image would retrieve the original input lattice.
// This requires that the input lattice have odd dimensions along
// each axis.
irtkRealImage * downsampleFactorTwo(irtkRealImage *in)
{

  double xspacing, yspacing, zspacing;
  double xaxis[3], yaxis[3], zaxis[3];

  irtkPoint origin;
  int xdim, ydim, zdim;

  in->GetPixelSize(&xspacing, &yspacing, &zspacing);
  origin = in->GetOrigin();
  in->GetOrientation(xaxis, yaxis, zaxis);

  xdim = in->GetX();
  ydim = in->GetY();
  zdim = in->GetZ();

  if (xdim % 2 == 0 || ydim % 2 == 0 || zdim % 2 == 0){
    cerr << "downsampleFactorTwo: Cannot apply to image with an even dimension." << endl;
    exit(0);
  }

  xdim = (xdim + 1) / 2;
  ydim = (ydim + 1) / 2;
  zdim = (zdim + 1) / 2;

  xspacing *= 2.0;
  yspacing *= 2.0;
  zspacing *= 2.0;
  irtkImageAttributes attr;
  attr._x = xdim;
  attr._y = ydim;
  attr._z = zdim;

  attr._dx = xspacing;
  attr._dy = yspacing;
  attr._dz = zspacing;

  attr._xorigin = origin._x;
  attr._yorigin = origin._y;
  attr._zorigin = origin._z;

  for (int i = 0; i < 3; ++i){
    attr._xaxis[i] = xaxis[i];
    attr._yaxis[i] = yaxis[i];
    attr._zaxis[i] = zaxis[i];
  }
  irtkRealImage * out = new irtkRealImage(attr);
//  irtkRealImage *out    = new irtkRealImage(xdim, ydim, zdim,
//                              xspacing, yspacing, zspacing,
//                              origin, xaxis, yaxis, zaxis);

  //   // Testing:
//   int i, j, k;
//   for (k = 0; k < zdim; ++k){
//     for (j = 0; j < ydim; ++j){
//       for (i = 0; i < xdim; ++i){
//         out->Put(i, j, k, in->Get(2*i, 2*j, 2*k));
//       }
//     }
//   }

  return out;
}


// Subdivide image in treating its values as a set of B-spline
// coefficients.  The output image needs to have dimensions 2x-1,
// 2y-1, 2z-1 where x, y, z are the dimensions of the input
// image.
void subdivide(irtkRealImage *in, irtkRealImage *out)
{
  int i, j, k, i1, j1, k1, i2, j2, k2;
  int xdim, ydim, zdim;
  double val;

  // Weights for subdivision
  double w[2][3];
  w[1][0] = 0;
  w[1][1] = 1.0/2.0;
  w[1][2] = 1.0/2.0;
  w[0][0] = 1.0/8.0;
  w[0][1] = 6.0/8.0;
  w[0][2] = 1.0/8.0;

  xdim = in->GetX();
  ydim = in->GetY();
  zdim = in->GetZ();

  if (out->GetX() != 2*xdim-1 || out->GetY() != 2*ydim-1 || out->GetZ() != 2*zdim-1){
    cerr << "subdivide: Image dimensions don't match. (" << xdim << ", " << ydim << ", " << zdim << ") ";
    cerr << out->GetX() << ", " << out->GetY() << ", " << out->GetZ() << ") " << endl;
    exit(1);
  }

  for (k = 1; k < zdim-2; ++k){
    for (j = 1; j < ydim-2; ++j){
      for (i = 1; i < xdim-2; ++i){

        for (k1 = 0; k1 < 2; ++k1){
	  for (j1 = 0; j1 < 2; ++j1){
            for (i1 = 0; i1 < 2; ++i1){

              val = 0;

              for (k2 = 0; k2 < 3; ++k2){
                for (j2 = 0; j2 < 3; ++j2){
                  for (i2 = 0; i2 < 3; ++i2){
                    val += w[i1][i2] * w[j1][j2] * w[k1][k2] * in->Get(i+i2-1, j+j2-1, k+k2-1);
                  }
                }
              }
              out->Put(2*i+i1, 2*j+j1, 2*k+k1, val);
	    }
	  }
	}
      }
    }
  }

}

// Get the intensity at a set of world coordinates for a lattice
// of coefficients.
double getIntensity(irtkRealImage *coeffs, double xin, double yin, double zin)
{
  //  double *xdata;
  irtkRealPixel *ptr2coeffs;

  double s, t, u, vi, vii, value; // B_I, B_J, B_K
  int lineOffset, sliceOffset, j, k, l, m, n, S, T, U;
  int xdim, ydim, zdim;
  double x, y, z;

  x = xin;
  y = yin;
  z = zin;

  coeffs->WorldToImage(x, y, z);

  xdim = coeffs->GetX();
  ydim = coeffs->GetY();
  zdim = coeffs->GetZ();

  // Now calculate the real stuff
  l = (int)floor(x);
  m = (int)floor(y);
  n = (int)floor(z);

  s = x - l;
  t = y - m;
  u = z - n;

  // Check bounds
  if (l - 1 < 0 || m - 1 < 0 || n - 1 < 0 ||
      l + 3 > xdim || m + 3 > ydim || n + 3 > zdim){
    return 0;
  }

  // Calculate offset
  //  i = (_x + 8) * (_y + 4);
  sliceOffset = xdim * ydim - 4 * xdim;
  lineOffset = xdim - 4;

  value = 0;

  S = round(DBL_LUTSIZE*s);
  T = round(DBL_LUTSIZE*t);
  U = round(DBL_LUTSIZE*u);

  //  xdata = &(_xdata[n-1][m-1][l-1]);
  ptr2coeffs = coeffs->GetPointerToVoxels(l - 1, m - 1, n - 1);

  for (k = 0; k < 4; ++k){
    vi = 0;

    for (j = 0; j < 4; ++j){
      // Inner most loop unrolled starts here
      vii  = *ptr2coeffs * LookupTable[S][0]; // B_I, i = 0
      ptr2coeffs++;

      vii += *ptr2coeffs * LookupTable[S][1]; // B_I, i = 1
      ptr2coeffs++;

      vii += *ptr2coeffs * LookupTable[S][2]; // B_I, i = 2
      ptr2coeffs++;

      vii += *ptr2coeffs * LookupTable[S][3]; // B_I, i = 3
      ptr2coeffs++;
      // Inner most loop unrolled stops here

      vi += vii * LookupTable[T][j]; // B_J
      ptr2coeffs += lineOffset;
    }

    value += vi * LookupTable[U][k]; // B_K
    ptr2coeffs += sliceOffset;
  }

  return value;
}

void evaluateCoeffsOnImageLattice(irtkRealImage *img, irtkRealImage *coeffs)
{
  int i, j, k;
  double x, y, z;
  int xdim, ydim, zdim;

  xdim = img->GetX();
  ydim = img->GetY();
  zdim = img->GetZ();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        x = i;
        y = j;
        z = k;
        img->ImageToWorld(x, y, z);
        img->Put(i, j, k, getIntensity(coeffs, x, y, z));
      }
    }
  }

}

void evaluateCoeffsOnImageLattice(irtkGreyImage *img, irtkRealImage *coeffs)
{
  int i, j, k;
  double x, y, z;
  int xdim, ydim, zdim;

  xdim = img->GetX();
  ydim = img->GetY();
  zdim = img->GetZ();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        x = i;
        y = j;
        z = k;
        img->ImageToWorld(x, y, z);
        img->Put(i, j, k, (int) round(getIntensity(coeffs, x, y, z)));
      }
    }
  }

}

void estimateCoeffs(int res)
{
  int i, j, k, l, m, n, inputIndex;
  int I, J, K;
  int xi, yi, zi, xdim, ydim, zdim;
  int inputXdim, inputYdim, inputZdim;
  double xw, yw, zw;
  double s, t, u;
  double value, phi, b, norm, tmp;

  xdim = _coeffs[res]->GetX();
  ydim = _coeffs[res]->GetY();
  zdim = _coeffs[res]->GetZ();

  cout << "Coefficients at level " << res + 1 << endl;
  cout << "Dimensions : " << xdim << ", " << ydim << ", " << zdim << endl;

  irtkRealPixel lo, hi;
  _coeffs[res]->GetMinMax(&lo, &hi);
  cout << "Min and max: " << lo << " " << hi << endl;


  // Loop over input images or packages.
  for (inputIndex = 0; inputIndex < _inputCount; ++inputIndex){

    cout << "Using image " << _input_name[inputIndex];
    cout << " with " << _dof_name[inputIndex] << " ... ";
    cout.flush();

    inputXdim = _input[inputIndex]->GetX();
    inputYdim = _input[inputIndex]->GetY();
    inputZdim = _input[inputIndex]->GetZ();

    // Loop over current input packages voxels.
    for (zi = 0; zi < inputZdim; ++zi){
      for (yi = 0; yi < inputYdim; ++yi){
        for (xi = 0; xi < inputXdim; ++xi){

          xw = xi;
          yw = yi;
          zw = zi;

          _input[inputIndex]->ImageToWorld(xw, yw, zw);
          _transf[inputIndex]->Transform(xw, yw, zw);

          value = _input[inputIndex]->Get(xi, yi, zi);
          value -= getIntensity(_coeffs[res], xw, yw, zw);

          _coeffs[res]->WorldToImage(xw, yw, zw);

          l = (int)floor(xw);
          m = (int)floor(yw);
          n = (int)floor(zw);

          if (l-1 < 0 || m-1 < 0 || n-1 < 0)
            continue;

          if (l+3 > xdim || m+3 > ydim || n+3 > zdim)
            continue;

          s = xw - l;
          t = yw - m;
          u = zw - n;

          norm = 0;
          for (k = 0; k < 4; ++k){
            for (j = 0; j < 4; ++j){
              for (i = 0; i < 4; ++i){
                tmp   = Bsp(i, s) * Bsp(j, t) * Bsp(k, u);
                norm += tmp * tmp;
              }
            }
          }

          for (k = 0; k < 4; ++k){
            K = n + k - 1;
            for (j = 0; j < 4; ++j){
              J = m + j - 1;
              for (i = 0; i < 4; ++i){
                I = l + i - 1;

                b = Bsp(i, s) * Bsp(j, t) * Bsp(k, u);
                phi = value * b / norm;

                _dx[res]->Put(I, J, K,
                                _dx[res]->Get(I, J, K)  + b * b * phi);

                _weights[res]->Put(I, J, K,
                                     _weights[res]->Get(I, J, K) + b * b);

              }
            }
          }

        }
      }
    }
    cout << " done." << endl;
  }

  for (zi = 0; zi < zdim; ++zi){
    for (yi = 0; yi < ydim; ++yi){
      for (xi = 0; xi < xdim; ++xi){

        if (_weights[res]->Get(xi, yi, zi) > 0){
          value = _coeffs[res]->Get(xi, yi, zi) +
            (_dx[res]->Get(xi, yi, zi) / _weights[res]->Get(xi, yi, zi));

          _coeffs[res]->Put(xi, yi, zi, value);
        }

      }
    }
  }

}


void usage()
{
  cerr << "Usage: multislicecorrect [coefficientLattice] [outputPrefix] [N] ";
  cerr << "[pack_1] .. [pack_N] [dof_1] .. [dof_N] <options>\n" << endl;
  cerr << endl;

  cerr << "\t[coefficientLattice] Image that defines the lattice on which the" << endl;
  cerr << "\t                     B-spline coefficients are estimated." << endl;
  cerr << "\t[outputPrefix]       Prefix to be used on the final reconstruction." << endl;
  cerr << "\t                     A '.nii.gz' suffix will be applied." << endl;
  cerr << "\t[pack_1] .. [pack_N] The input packages." << endl;
  cerr << "\t[dof_1]  .. [dof_N]  The transformations between the input packages" << endl;
  cerr << "\t                     and the coefficient lattice.  In each case, the" << endl;
  cerr << "\t                     package should be the target of the transfomation." << endl;
  cerr << "\t                     The string \"id\" can be used for any of [dof_1] ... [dof_N]" << endl;
  cerr << "\t                     In this case an identity transformation is used." << endl;
  cerr << "\t" << endl;
  cerr << "Options:" << endl;
  cerr << "\t-levels [number] :   Number of resolution levels required (Default 3)." << endl;
  cerr << "\t                     These are estimated in a coarse to fine fashion " << endl;
  cerr << "\t                     in reverse order, e.g. for 3 levels:" << endl;
  cerr << "\t                     Coarsest level 3 is estimated first followed by" << endl;
  cerr << "\t                     finer level 2 then finest resolution level 1." << endl;
  cerr << "\t-lastlevel [n]       The final level to be estimated." << endl;
  cerr << "\t-reconLattice [img]  Evaluate the reconstruction estimate on the lattice of " << endl;
  cerr << "\t                     image [img].  The default is to use the lattice" << endl;
  cerr << "\t                     of image [coefficientLattice]." << endl;
  cerr << "\t-allLevels           Write the reconstruction output for all estimated" << endl;
  cerr << "\t                     levels.  The output images will show a number to\t" << endl;
  cerr << "\t                     indicate the level of reconstruction." << endl;
  cerr << "\t" << endl;
  cerr << "\t" << endl;

  cerr << "\t" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, inputIndex, ok, pad;
  char buffer[256];

  // Default values.
  int levels = 3;
  int lastLevel = 1;
  int writeCoeffs = False;
  int writeWeights = False;
  int writeAllLevels = False;

  if (argc < 6)
    usage();

  _template_name = argv[1];
  argc--;
  argv++;

  _output_name = argv[1];
  argc--;
  argv++;

  _inputCount = atoi(argv[1]);
  argc--;
  argv++;

  // Parse input file and transformation names
  _input_name = new char *[MAX_IMAGES];
  _dof_name   = new char *[MAX_IMAGES];

  for (inputIndex = 0; inputIndex < _inputCount; ++inputIndex){
    _input_name[inputIndex] =  argv[1];
    argc--;
    argv++;
  }

  for (inputIndex = 0; inputIndex < _inputCount; ++inputIndex){
    _dof_name[inputIndex] =  argv[1];
    argc--;
    argv++;
  }

  // Parse options.
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-weights") == 0)){
      argc--;
      argv++;
      writeWeights = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-coeffs") == 0)){
      argc--;
      argv++;
      writeCoeffs = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-levels") == 0)){
      argc--;
      argv++;
      levels = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-lastlevel") == 0)){
      argc--;
      argv++;
      lastLevel = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-reconLattice") == 0)){
      argc--;
      argv++;
      _evalTarget_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-allLevels") == 0)){
      argc--;
      argv++;
      writeAllLevels = True;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (lastLevel < 1 || lastLevel > levels){
    cerr << "Last level outside possible range given by number of levels." << endl;
    exit(1);
  }

  cout << "Downsampling to level : " << levels << endl;
  cout << "Last level to be optimised : " << lastLevel << endl;


  // Read images and dofs:
  _input = new irtkRealImage*[_inputCount];
  _transf = new irtkTransformation*[_inputCount];

  for (i = 0; i < _inputCount; ++i){
    cout << "Reading image " << _input_name[i] << " ... " ;
    cout.flush();
    _input[i] = new irtkRealImage(_input_name[i]);

    if (strcmp(_dof_name[i], "id") == 0){
      _transf[i] = new irtkRigidTransformation;
      cout << "using identity dof." << endl;
    } else {
      cout << "using input dof " << _dof_name[i] << " ... "; cout.flush();
      _transf[i] = irtkTransformation::New(_dof_name[i]);
      cout << "done" << endl;
    }
  }

  // Final resolution template.
  _template  = new irtkRealImage(_template_name);

  if (_evalTarget_name != NULL){
    _evalTarget = new irtkRealImage(_evalTarget_name);
  } else {
    _evalTarget = new irtkRealImage(_template_name);
  }

  // Initialisation:
  initialiseLookupTable();
  _coeffs  = new irtkRealImage*[levels];
  _dx      = new irtkRealImage*[levels];
  _weights = new irtkRealImage*[levels];


  // Put zeros around the image, enough for the subdivsions
  // required.
  pad = 1;
  for (i = 0; i < levels; ++i){
    pad *= 2;
  }
  _coeffs[0] = zeroPad(_template, pad);

  // Make the dimensions right for subdivision.
  _temp = growImage(_coeffs[0], levels);
  delete _coeffs[0];
  _coeffs[0] = _temp;


  // Make the pyramid.
  for (i = 1; i < levels; ++i){
    _coeffs[i] = downsampleFactorTwo(_coeffs[i - 1]);
  }

  // Grab memory for _dx and _weights
  for (i = 0; i < levels; ++i){
    _dx[i]      = new irtkRealImage(*_coeffs[i]);
    _weights[i] = new irtkRealImage(*_coeffs[i]);
    clearRealImage(_coeffs[i]);
    clearRealImage(_dx[i]);
    clearRealImage(_weights[i]);
  }

  // Estimate coefficients, starting with coarsest lattice.
  for (i = levels - 1; i >= lastLevel - 1; --i){

    cout << "Estimating coefficients for level " << i+1 << ":" << endl;
    estimateCoeffs(i);
    cout << "Level " << i+1 << " done" << endl;
    cout << "-----" << endl;

    // Real
//     _temp = new irtkRealImage(*_evalTarget);
//     clearRealImage(_temp);
//     evaluateCoeffsOnImageLattice(_temp, _coeffs[i]);
//     sprintf(buffer, "%s-level.%d.nii.gz", _output_name, i+1);
//     _temp->Write(buffer);
//     delete _temp;

    // Grey
    _temp2 = new irtkGreyImage(*_evalTarget);
    clearGreyImage(_temp2);
    evaluateCoeffsOnImageLattice(_temp2, _coeffs[i]);

    if (i == lastLevel - 1 && writeAllLevels == False){
      sprintf(buffer, "%s.nii.gz", _output_name);
      _temp2->Write(buffer);
    } else if (writeAllLevels == True){
      sprintf(buffer, "%s-level.%d.nii.gz", _output_name, i+1);
      _temp2->Write(buffer);
    }

    delete _temp2;

    // Subdivide for the next size lattice if needed.
    if (i > lastLevel - 1){
      cout << "Subdividing ... ";
      cout.flush();
      subdivide(_coeffs[i], _coeffs[i - 1]);
      cout << "done." << endl;
    }

    if (writeWeights == True){
      sprintf(buffer, "%s-weights_%d.nii.gz", _output_name, i+1);
      _weights[i]->Write(buffer);
    }

    if (writeCoeffs == True){
      sprintf(buffer, "%s-coeffs_%d.nii.gz", _output_name, i+1);
      _coeffs[i]->Write(buffer);
    }

  }

}



//////////////////////////////////////////////////////////////////
// Graveyard.
//////////////////////////////////////////////////////////////////

//   double k1 = 2.0/3.0;
//   double k2 = 1.0/6.0;
//   double bspReduced[3];
//   bspReduced[0] = bspReduced[2] = k2;
//   bspReduced[1] = k1;


void conv_3D(irtkRealImage *image, irtkRealImage *temp, double *ker, int lKern)
{
  int i, j, k, t;
  int xdim, ydim, zdim;
  int xstart, xend, ystart, yend, zstart, zend;
  int midKern;
  double sum;
  irtkRealPixel *imgPtr, *tempPtr, *offsetPtr;
  int xstep, ystep, zstep;

  // 1. lKern should be odd.
  // 2. image and temp must have same dimensions.
  // Should check.
  // temp ends up being the output file.
  // the image 'image' is affected by this convolution :-(

  xdim = image->GetX();
  ydim = image->GetY();
  zdim = image->GetZ();

  midKern = -1 + (1 + lKern) / 2;

  xstart = ystart = zstart = midKern;
  xend = xdim - 1 - midKern;
  yend = ydim - 1 - midKern;
  zend = zdim - 1 - midKern;

  xstep = 1;
  ystep = xdim;
  zstep = xdim * ydim;

  // Convolve along x.

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      imgPtr  = image->GetPointerToVoxels(xstart, j, k);
      tempPtr = temp->GetPointerToVoxels(xstart, j, k);

      for (i = xstart; i <= xend; ++i){
        sum = 0.0;
        for (t = 0; t < lKern; ++t){
          offsetPtr = imgPtr + (t - midKern) * xstep;
          sum += ker[t] * (*offsetPtr);
        }
        *tempPtr = sum;
        tempPtr = tempPtr + xstep;
        imgPtr  = imgPtr  + xstep;
      }
    }
  }

  // Convolve along y. Swap roles of image and temp.
  for (k = 0; k < zdim; ++k){
    for (i = 0; i < xdim; ++i){
      imgPtr  = image->GetPointerToVoxels(i, ystart, k);
      tempPtr = temp->GetPointerToVoxels(i, ystart, k);

      for (j = ystart; j <= yend; ++j){
        sum = 0.0;
        for (t = 0; t < lKern; ++t){
          offsetPtr = tempPtr + (t - midKern) * ystep;
          sum += ker[t] * (*offsetPtr);
        }
        *imgPtr = sum;
        tempPtr = tempPtr + ystep;
        imgPtr  = imgPtr  + ystep;
      }
    }
  }

  // Convolve along z.
  for (j = 0; j < ydim; ++j){
    for (i = 0; i < xdim; ++i){
      imgPtr  = image->GetPointerToVoxels(i, j, zstart);
      tempPtr = temp->GetPointerToVoxels(i, j, zstart);

      for (k = zstart; k <= zend; ++k){
        sum = 0.0;
        for (t = 0; t < lKern; ++t){
          offsetPtr = imgPtr + (t - midKern) * zstep;
          sum += ker[t] * (*offsetPtr);
        }
        *tempPtr = sum;
        tempPtr = tempPtr + zstep;
        imgPtr  = imgPtr  + zstep;
      }
    }
  }
}

void conv_3D_long(irtkRealImage *image, irtkRealImage *temp, double *ker, int lKern)
{
  int i, j, k, s, t, u;
  int xdim, ydim, zdim;
  int xstart, xend, ystart, yend, zstart, zend;
  int midKern;
  double sum;

  // 1. lKern should be odd.
  // 2. image and temp must have same dimensions.
  // Should check.
  // temp ends up being the output file.
  // the image 'image' is affected by this convolution :-(

  xdim = image->GetX();
  ydim = image->GetY();
  zdim = image->GetZ();

  midKern = -1 + (1 + lKern) / 2;

  xstart = ystart = zstart = midKern;
  xend = xdim - 1 - midKern;
  yend = ydim - 1 - midKern;
  zend = zdim - 1 - midKern;

  for (k = zstart; k <= zend; ++k){
    for (j = ystart; j <= yend; ++j){
      for (i = xstart; i <= xend; ++i){
        sum = 0.0;
        for (u = 0; u < lKern; ++u){
          for (t = 0; t < lKern; ++t){
            for (s = 0; s < lKern; ++s){
              sum += ker[u]*ker[t]*ker[s]*(image->Get(i + s - midKern, j + t - midKern, k + u - midKern));
//               if (i == 10 && j == 10 && k == 7){
//                 cout << i + s - midKern << " " <<  j + t - midKern << " " <<  k + u - midKern << " : ";
//                 cout << image->Get(i + s - midKern, j + t - midKern, k + u - midKern) << " : ";
//                 cout << ker[s] << " " << ker[t] << " " << ker[u] << endl;
//               }
            }
          }
        }
        temp->Put(i, j, k, sum);
      }
    }
  }


}



double getIntensity_verbose(irtkRealImage *coeffs, double x, double y, double z)
{
  //  double *xdata;
  irtkRealPixel *ptr2coeffs;

  double s, t, u, B_J, B_K, vi, vii, value; // B_I
  int lineOffset, sliceOffset, j, k, l, m, n, S, T, U;
  int xdim, ydim, zdim;

  coeffs->WorldToImage(x, y, z);

  xdim = coeffs->GetX();
  ydim = coeffs->GetY();
  zdim = coeffs->GetZ();

  // Now calculate the real stuff
  l = (int)floor(x);
  m = (int)floor(y);
  n = (int)floor(z);
  s = x-l;
  t = y-m;
  u = z-n;

  // Check bounds
  if (l - 1 < 0 || m - 1 < 0 || n - 1 < 0 ||
      l + 3 > xdim || m + 3 > ydim || n + 3 > zdim){
    return 0;
  }


  // Calculate offset
  //  i = (_x + 8) * (_y + 4);
  sliceOffset = xdim * ydim - 4 * xdim;
  lineOffset = xdim - 4;

  value = 0;

  S = round(DBL_LUTSIZE*s);
  T = round(DBL_LUTSIZE*t);
  U = round(DBL_LUTSIZE*u);

  //  xdata = &(_xdata[n-1][m-1][l-1]);
  ptr2coeffs = coeffs->GetPointerToVoxels(l - 1, m - 1, n - 1);

  for (k = 0; k < 4; ++k){
    B_K = LookupTable[U][k];
    vi = 0;

    for (j = 0; j < 4; ++j){
      B_J = LookupTable[T][j];

      // Inner most loop unrolled starts here
      vii  = *ptr2coeffs * LookupTable[S][0];
      ++ptr2coeffs;

//       cout << *ptr2coeffs << " : " << LookupTable[S][0] << " " << B_J << " " << B_K << endl;

      vii += *ptr2coeffs * LookupTable[S][1];
      ++ptr2coeffs;

//       cout << *ptr2coeffs << " : " << LookupTable[S][1] << " " << B_J << " " << B_K << endl;

      vii += *ptr2coeffs * LookupTable[S][2];
      ++ptr2coeffs;

//       cout << *ptr2coeffs << " : " << LookupTable[S][2] << " " << B_J << " " << B_K << endl;

      vii += *ptr2coeffs * LookupTable[S][3];
      ++ptr2coeffs;

//       cout << *ptr2coeffs << " : " << LookupTable[S][3] << " " << B_J << " " << B_K << endl;

      // Inner most loop unrolled stops here

      vi += vii * B_J;
      ptr2coeffs += lineOffset;
    }

    value += vi * B_K;
    ptr2coeffs += sliceOffset;

  }

  ////////////////////////////////////////////////////
  int i, I, J, K;
  double tmp = 0;
  for (k = 0; k < 4; ++k){
    K = n - 1 + k;
    for (j = 0; j < 4; ++j){
      J = m - 1 + j;
      for (i = 0; i < 4; ++i){
        I = l - 1 + i;
        tmp += LookupTable[S][i] * LookupTable[T][j] * LookupTable[U][k] * coeffs->Get(I, J, K);
      }
    }
  }
  cout << value << " \t" << tmp << endl;
  /////////////////////////////////////////////////////


  return value;

}

// Upsample a lattice in the FFD style.
// Downsample in way such that FFD style subdivision applied to
// the output image would retrieve the original input lattice.
// This requires that the input lattice have odd dimensions along
// each axis.
irtkRealImage * upsampleFactorTwo(irtkRealImage *in){

  double xspacing, yspacing, zspacing;
  double xaxis[3], yaxis[3], zaxis[3];

  irtkPoint origin;
  int xdim, ydim, zdim;

  in->GetPixelSize(&xspacing, &yspacing, &zspacing);

  origin = in->GetOrigin();
  in->GetOrientation(xaxis, yaxis, zaxis);

  xdim = in->GetX();
  ydim = in->GetY();
  zdim = in->GetZ();

  xdim = 2 * xdim - 1;
  ydim = 2 * ydim - 1;
  zdim = 2 * zdim - 1;

  xspacing /= 2.0;
  yspacing /= 2.0;
  zspacing /= 2.0;

  irtkImageAttributes attr;
  attr._x = xdim;
  attr._y = ydim;
  attr._z = zdim;

  attr._dx = xspacing;
  attr._dy = yspacing;
  attr._dz = zspacing;

  attr._xorigin = origin._x;
  attr._yorigin = origin._y;
  attr._zorigin = origin._z;

  for (int i = 0; i < 3; ++i){
    attr._xaxis[i] = xaxis[i];
    attr._yaxis[i] = yaxis[i];
    attr._zaxis[i] = zaxis[i];
  }
  irtkRealImage *out = new irtkRealImage(attr);

//  irtkRealImage *out    = new irtkRealImage(xdim, ydim, zdim,
//                              xspacing, yspacing, zspacing,
//                              origin, xaxis, yaxis, zaxis);

//   // Testing:
//   int i, j, k;
//   xdim = in->GetX();
//   ydim = in->GetY();
//   zdim = in->GetZ();
//   for (k = 0; k < zdim; ++k){
//     for (j = 0; j < ydim; ++j){
//       for (i = 0; i < xdim; ++i){
//         out->Put(2*i, 2*j, 2*k, in->Get(i, j, k));
//       }
//     }
//   }

  return out;

}


