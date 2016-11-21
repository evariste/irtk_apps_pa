
#include <irtkImage.h>

// Code that strips down the IRTK B-spline image interpolation
// to 1-D so that it can be applied to a single array.
//
// See main function where it is applied to some toy example data.
//
// Any bugs etc. please let me know, thanks, Paul.

// Coefficient helper functions.
double InitialCausalCoefficient(double c[], int DataLength, double z, double Tolerance);
double InitialAntiCausalCoefficient(double c[], int DataLength, double z);
void ConvertToInterpolationCoefficients(double *c, int DataLength, double *z, int NbPoles, double Tolerance);
void ComputeCoefficients(double *data, int _SplineDegree, int nSamples);

// Evaluate at a specific grid location.
double Evaluate(double u_x, double *_coeff, int _SplineDegree, int nSamples);

// Interpolate some data.
void bspline_interp_1d_pa(int nSamples, double *ySamp, int splineDegree,
                    int nRequired, double *u_x, double *y_interp);
// nSamples: How many samples in the measured data
// ySamp: The measurements themselves.
// splineDegree: Can be from 2 to 5
// nRequired: How many points do we want to interpolate at?
// u_x: The grid coordinates where we want to interpolate at
// y_interp: Where to store the interpolated values.



////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  int i;

  // The x range where the samples are
  double xMin = 0;
  double xMax = 5;

  // How many points to sample (including end points)
  int nSamples = 11;

  // step size for the sampled points
  double h = (xMax - xMin) / (nSamples - 1);

  // Make up some sample data.
  double *xSamp, *ySamp;
  xSamp = new double[nSamples];
  ySamp = new double[nSamples];

  for (i = 0; i < nSamples; ++i){
    double xi = xMin + h * i;
    xSamp[i] = xi;

    // Some function to generate y data (taken from GSL example).
    ySamp[i] = cos(xi) * exp(-0.1 * xi);
  }


  // Print the data
  for (i = 0; i < nSamples; ++i){
    cout << xSamp[i] << "," << ySamp[i] << endl;
  }
  cout << "-------------------------" << endl;


  // The step size for the points where we want to do a finer
  // grained interpolation.
  double hSmall = 0.25;

  // How many interpolated points will there be with the finer
  // step size?
  int nRequired = round((xMax - xMin) / hSmall);

  if ((xMin + hSmall * nRequired) > xMax)
    nRequired -= 1;

  // Storage for locations where we require to interpolate -
  // and the corresponding estimates.
  double *x_req, *y_interp;
  x_req = new double[nRequired];
  y_interp = new double[nRequired];
  // Locations:
  for (i = 0; i < nRequired; i++){
    double xi = xMin + i * hSmall;
    x_req[i] = xi;
  }

  // Grid coordinates for the locations where we want to
  // interpolate - using the original sampling step, h.
  double *u_x = new double[nRequired];
  for (i = 0; i < nRequired; i++){
    u_x[i] = (x_req[i] - xMin) / h;
  }

  // Now actually do the interpolation.
  int splineDegree = 3;
  bspline_interp_1d_pa(nSamples, ySamp, splineDegree, nRequired, u_x, y_interp);


  // Print the interpolated values
  for (i = 0; i < nRequired; i++){
    cout << x_req[i] << "," << y_interp[i] << endl;
  }
  cout << "-------------------------" << endl;

  // Clean up.
  delete [] xSamp;
  delete [] ySamp;
  delete [] x_req;
  delete [] y_interp;
  delete [] u_x;

}

////////////////////////////////////////////////////////

void bspline_interp_1d_pa(int nSamples, double *ySamp, int splineDegree,
                    int nRequired, double *u_x, double *y_interp)
{
  // nSamples: How many samples in the measured data
  // ySamp: The measurements themselves.
  // splineDegree: Can be from 2 to 5
  // nRequired: How many points do we want to interpolate at?
  // u_x: The grid coordinates where we want to interpolate at
  // y_interp: Where to store the interpolated values.

  int i;
  double *coeff = new double[nSamples];

  for (i = 0; i < nSamples; ++i){
    coeff[i] = ySamp[i];
  }

  ComputeCoefficients(coeff, splineDegree, nSamples);

  for (i = 0; i < nRequired; i++){
    y_interp[i] = Evaluate(u_x[i], coeff, splineDegree, nSamples);
  }

  delete [] coeff;

}

////////////////////////////////////////////////////////

double Evaluate(double u_x, double *_coeff, int _SplineDegree, int nSamples)
{
  int i, m;
  int index[6];
  double weight[6];
  double value, w, w2, w4, t, t0, t1;

  int nHalf = 2 * nSamples - 2;

  /* compute the interpolation indexes */
  if (_SplineDegree & 1) {
    i = (int)floor(u_x) - _SplineDegree / 2;

    for (m = 0; m <= _SplineDegree; m++) {
      index[m] = i++;
    }
  } else {
    i = (int)floor(u_x + 0.5) - _SplineDegree / 2;

    for (m = 0; m <= _SplineDegree; m++) {
      index[m] = i++;
    }
  }

  /* compute the interpolation weights */
  switch (_SplineDegree) {
  case 2:
    /* u_x */
    w = u_x - (double)index[1];
    weight[1] = 3.0 / 4.0 - w * w;
    weight[2] = (1.0 / 2.0) * (w - weight[1] + 1.0);
    weight[0] = 1.0 - weight[1] - weight[2];
    break;
  case 3:
    /* u_x */
    w = u_x - (double)index[1];
    weight[3] = (1.0 / 6.0) * w * w * w;
    weight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - weight[3];
    weight[2] = w + weight[0] - 2.0 * weight[3];
    weight[1] = 1.0 - weight[0] - weight[2] - weight[3];
    break;
  case 4:
    /* u_x */
    w = u_x - (double)index[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    weight[0] = 1.0 / 2.0 - w;
    weight[0] *= weight[0];
    weight[0] *= (1.0 / 24.0) * weight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    weight[1] = t1 + t0;
    weight[3] = t1 - t0;
    weight[4] = weight[0] + t0 + (1.0 / 2.0) * w;
    weight[2] = 1.0 - weight[0] - weight[1] - weight[3] - weight[4];
    break;
  case 5:
    /* u_x */
    w = u_x - (double)index[2];
    w2 = w * w;
    weight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    weight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - weight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    weight[2] = t0 + t1;
    weight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    weight[1] = t0 + t1;
    weight[4] = t0 - t1;
    break;
  default:
    printf("Invalid spline degree\n");
    return(0.0);
  }

  /* apply the mirror boundary conditions */
  for (m = 0; m <= _SplineDegree; m++) {
    index[m] = (nSamples == 1) ? (0) : ((index[m] < 0) ? (-index[m] - nHalf * ((-index[m]) / nHalf)) : (index[m] - nHalf * (index[m] / nHalf)));
    if (nSamples <= index[m]) {
      index[m] = nHalf - index[m];
    }
  }

  /* perform interpolation */

  value = 0.0;

  for (i = 0; i <= _SplineDegree; i++) {
    value += weight[i] * _coeff[index[i]];
  }

  return(value);
}

////////////////////////////////////////////////////////

double InitialCausalCoefficient(double c[], int DataLength, double z, double Tolerance)
{
  double Sum, zn, z2n, iz;
  int n, Horizon;

  /* this initialization corresponds to mirror boundaries */
  Horizon = DataLength;
  if (Tolerance > 0.0) {
    Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
  }
  if (Horizon < DataLength) {
    /* accelerated loop */
    zn = z;
    Sum = c[0];
    for (n = 1; n < Horizon; n++) {
      Sum += zn * c[n];
      zn *= z;
    }
    return(Sum);
  } else {
    /* full loop */
    zn = z;
    iz = 1.0 / z;
    z2n = pow(z, (double)(DataLength - 1));
    Sum = c[0] + z2n * c[DataLength - 1];
    z2n *= z2n * iz;
    for (n = 1; n <= DataLength - 2; n++) {
      Sum += (zn + z2n) * c[n];
      zn *= z;
      z2n *= iz;
    }
    return(Sum / (1.0 - zn * zn));
  }
}

////////////////////////////////////////////////////////


double InitialAntiCausalCoefficient(double c[], int DataLength, double z)
{
  /* this initialization corresponds to mirror boundaries */
  return((z / (z * z - 1.0)) * (z * c[DataLength - 2] + c[DataLength - 1]));
}

////////////////////////////////////////////////////////

void ConvertToInterpolationCoefficients(double *c, int DataLength, double *z, int NbPoles, double Tolerance)
{
  double Lambda = 1.0;
  int n, k;

  /* special case required by mirror boundaries */
  if (DataLength == 1) {
    return;
  }

  /* compute the overall gain */
  for (k = 0; k < NbPoles; k++) {
    Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
  }

  /* apply the gain */
  for (n = 0; n < DataLength; n++) {
    c[n] *= Lambda;
  }

  /* loop over all poles */
  for (k = 0; k < NbPoles; k++) {
    /* causal initialization */
    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
    /* causal recursion */
    for (n = 1; n < DataLength; n++) {
      c[n] += z[k] * c[n - 1];
    }
    /* anticausal initialization */
    c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
    /* anticausal recursion */
    for (n = DataLength - 2; 0 <= n; n--) {
      c[n] = z[k] * (c[n + 1] - c[n]);
    }
  }
}

////////////////////////////////////////////////////////

void ComputeCoefficients(double *data, int _SplineDegree, int nSamples)
{
//  double *data;
  double Pole[2];
  int NbPoles;

  /* recover the poles from a lookup table */
  switch (_SplineDegree) {
  case 2:
    NbPoles = 1;
    Pole[0] = sqrt(8.0) - 3.0;
    break;
  case 3:
    NbPoles = 1;
    Pole[0] = sqrt(3.0) - 2.0;
    break;
  case 4:
    NbPoles = 2;
    Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
    Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
    break;
  case 5:
    NbPoles = 2;
    Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
              - 13.0 / 2.0;
    Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
              - 13.0 / 2.0;
    break;
  default:
    cerr << "Invalid spline degree" << endl;
    exit(1);
  }

  /* convert the image samples into interpolation coefficients */
  ConvertToInterpolationCoefficients(data, nSamples, Pole, NbPoles, DBL_EPSILON);

}





