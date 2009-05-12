#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
// Estimate the plane of symmetry for a mesh (e.g. an extracted cortex).

#include <irtkImage.h>
#include <irtkLocator.h>

#include <nr.h>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkMath.h>
#include <vtkPlaneSource.h>

#ifndef M_PI_2
# define M_PI_2   1.57079632679489661923  // pi/2
#endif


float symmetryMeasure(float *normal, float *centre);
void normal2phiTheta(float *normal, float &phi, float &theta);
void phiTheta2normal(float phi, float theta, float* normal);
float symmetryMeasure(float phi, float theta, float *centre);
float objectiveFuncForNR(float params[]);

char *surface_name = NULL;

vtkPolyData* _surface;
irtkLocator *_locator;

int main(int argc, char **argv){

  int ok;
  int noOfPoints;
  double cx, cy, cz;
  double pt[3];
  float normal[3], centre[3];
  double xx, xy, xz, yy, yz, zz;
  int i, j;
  float phi, theta;

  double val;

  if (argc < 2){
    cerr << argv[0] << " [surface] " << endl;
    cerr << "" << endl;
    exit(1);
  }

  surface_name = argv[1];
  argv++;
  argc--;


  while (argc > 1){
    ok = False;
//     if ((ok == False) && (strcmp(argv[1], "-option") == 0)){
//       argc--;
//       argv++;
//       ok = True;
//     }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(surface_name);
  reader->Update();

  _surface = vtkPolyData::New();
  _surface = reader->GetOutput();

  noOfPoints = _surface->GetNumberOfPoints();

  if (noOfPoints <= 0){
    cerr << "Surface has no points." << endl;
    exit(1);
  }

  // Find centre of gravity.
  cx = cy = cz = 0.0;
  for (int i = 0; i < noOfPoints; i++){
    _surface->GetPoint (i, pt);
    cx += pt[0];
    cy += pt[1];
    cz += pt[2];
  }

  cx /= noOfPoints;
  cy /= noOfPoints;
  cz /= noOfPoints;

  cout << "Centre of gravity : (" << cx << ", " << cy << ", " << cz << ")" <<  endl;

  xx = xy = xz = yy = yz = zz = 0.0;

  irtkMatrix cov;
  cov.Initialize(3, 3);


  for (int i = 0; i < noOfPoints; i++){
    _surface->GetPoint (i, pt);
    xx += (pt[0] - cx) * (pt[0] - cx);
    xy += (pt[0] - cx) * (pt[1] - cy);
    xz += (pt[0] - cx) * (pt[2] - cz);
    yy += (pt[1] - cy) * (pt[1] - cy);
    yz += (pt[1] - cy) * (pt[2] - cz);
    zz += (pt[2] - cz) * (pt[2] - cz);
  }

  xx /= noOfPoints;
  xy /= noOfPoints;
  xz /= noOfPoints;
  yy /= noOfPoints;
  yz /= noOfPoints;
  zz /= noOfPoints;

  cout << "Finding covariance" << endl;

  cov(0, 0) = xx;
  cov(1, 1) = yy;
  cov(2, 2) = zz;

  cov(0, 1) = xy;
  cov(1, 0) = xy;

  cov(0, 2) = xz;
  cov(2, 0) = xz;

  cov(1, 2) = yz;
  cov(2, 1) = yz;

  cout << "Finding eigenstuff " << endl;

  irtkMatrix evecs;
  irtkVector evals;


  cov.Eigenvalues(evecs, evals);

  cout << "Covariance:" << endl;
  cov.Print();
  cout << " " << endl;

  cout << "Evecs: " << endl;
  evecs.Print();
  cout << " " << endl;

  cout << "Evals: " << endl;
  evals.Print();
  cout << " " << endl;

  // Set up locator.
  _locator = new irtkLocator;
  _locator->SelectLocatorType(1);
  _locator->SetDataSet(_surface);

  centre[0] = cx;
  centre[1] = cy;
  centre[2] = cz;

  float measure, minMeasure;
  int minCol = 0;

  float params[6];

  minMeasure = FLT_MAX;

  for (j = 0; j < 3; ++j){
    for (i = 0; i < 3; ++i){
      normal[i] = evecs(i, j);
    }
    cout << "Symmetry measure for e-vec " << j + 1 << ": ";
    measure = symmetryMeasure(normal, centre);
    cout << measure << endl;
    if (measure < minMeasure){
      minMeasure = measure;
      minCol = j;
    }
  }

  // Assign minimising e-vec to the normal.
  for (i = 0; i < 3; ++i){
    normal[i] = evecs(i, minCol);
  }

  cout << "Eigenvector " << minCol + 1 << " gives minimum symmetry measure." << endl;
  cout << "Normal      : " << normal[0] << " " << normal[1] << " " << normal[2] << endl;

  // Express normal in degrees:
  normal2phiTheta(normal, phi, theta);

  // Set up parameters for NR call.
  params[1] = phi;
  params[2] = theta;
  params[3] = centre[0];
  params[4] = centre[1];
  params[5] = centre[2];

  // Initial set of directions (unit vectors);
  float **initDirs;
  initDirs = new float* [6];

  for (i = 0; i < 6; ++i){
    initDirs[i] = new float[6];

    for (j = 0; j < 6; ++j){
      initDirs[i][j] = 0.0;
    }
    initDirs[i][i] = 1.0;
  }

  // Tolerance for stopping condition.
  float ftol;
  ftol = 0.5;

  // Number of iterations taken.
  int iters;

  // Final value of objective function.
  float fReturned;

  cout << "Powell minimisation ...";
  cout.flush();
  powell(params, initDirs, 5, ftol, &iters, &fReturned, objectiveFuncForNR);
  cout << "... done." << endl;

  for (i = 1; i < 6; ++i){
    cout << params[i] << endl;
  }

  phi = params[1];
  theta = params[2];
  centre[0] = params[3];
  centre[1] = params[4];
  centre[2] = params[5];
  phiTheta2normal(phi, theta, normal);

  cout << "After optimisation: " << endl;
  cout << "Normal       : " << normal[0] << " " << normal[1] << " " << normal[2] << endl;
  cout << "Plane centre : " << centre[0] << " " << centre[1] << " " << centre[2] << endl;

  cout << "Measure      : " << objectiveFuncForNR(params) << endl;

  // Now make a plane.
  double surfaceBounds[6];
  double xmin, ymin, zmin, xmax, ymax, zmax;

  _surface->GetBounds(surfaceBounds);
  xmin = surfaceBounds[0];
  ymin = surfaceBounds[1];
  zmin = surfaceBounds[2];
  xmax = surfaceBounds[3];
  ymax = surfaceBounds[4];
  zmax = surfaceBounds[5];


  vtkPlaneSource *plane = vtkPlaneSource::New();
  plane->SetResolution(10, 10);
  // Find a couple of points in the plane.

  double nDotC;
  nDotC = vtkMath::Dot(normal, centre);
  // Try and get points for the plane that cover the object
  if (fabs(normal[2]) > 0.05){
    val = (nDotC - normal[0]*xmin - normal[1]*ymax) / normal[2];
    plane->SetOrigin(xmin, ymax, val);
    val = (nDotC - normal[0]*xmin - normal[1]*ymin) / normal[2];
    plane->SetPoint1(xmin, ymin, val);
    val = (nDotC - normal[0]*xmax - normal[1]*ymax) / normal[2];
    plane->SetPoint2(xmax, ymax, val);
  } else if (fabs(normal[1]) > 0.05){
    val = (nDotC - normal[0]*xmin - normal[2]*zmax) / normal[1];
    plane->SetOrigin(xmin, val, zmax);
    val = (nDotC - normal[0]*xmin - normal[2]*zmin) / normal[1];
    plane->SetPoint1(xmin, val, zmin);
    val = (nDotC - normal[0]*xmax - normal[2]*zmax) / normal[1];
    plane->SetPoint2(xmax, val, zmax);
  } else {
    val = (nDotC - normal[1]*ymin - normal[2]*zmax) / normal[0];
    plane->SetOrigin(val, ymin, zmax);
    val = (nDotC - normal[1]*ymin - normal[2]*zmin) / normal[0];
    plane->SetPoint1(val, ymin, zmin);
    val = (nDotC - normal[1]*ymax - normal[2]*zmax) / normal[0];
    plane->SetPoint2(val, ymax, zmax);
  }

  vtkPolyData *output = vtkPolyData::New();
  output = plane->GetOutput();
  output->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(output);
  writer->SetFileName("plane.vtk");
  writer->Update();
  writer->Write();



  for (i = 0; i < 6; ++i){
      delete [] initDirs[i];
  }
  delete [] initDirs;

}

float objectiveFuncForNR(float params[])
{
  //                     phi        theta      centre
  float phi, theta;
  float centre[3];
  float measure;

  // Pesky NR one-indexing.
  phi = params[1];
  theta = params[2];
  centre[0] = params[3];
  centre[1] = params[4];
  centre[2] = params[5];

  measure = symmetryMeasure(phi, theta, centre);

//  cout << "  " << phi << " " << theta << " " << centre[0] << " " << centre[1] << " " << centre [2] << " " << measure << endl;

  return measure;
}

float symmetryMeasure(float phi, float theta, float *centre)
{
  float normal[3];
  phiTheta2normal(phi, theta, normal);
  return symmetryMeasure(normal, centre);
}

// For each point in _surface, reflect in the plane defined by the arguments and find the distance
// to the surface.  The mean distance over all points is returned.
float symmetryMeasure(float *normal, float *centre)
{
  int i, j, noOfPoints;
  double pt[3], ptTemp[3];
  float diff[3];
  double temp, compNorm;
  double sumDist;

  // Ensure unit normal.
  temp = vtkMath::Norm(normal);
  for (j = 0; j < 3; ++j)
    normal[j] /= temp;

  noOfPoints = _surface->GetNumberOfPoints();
  sumDist = 0.0;

  for (i = 0; i < noOfPoints; ++i){
    _surface->GetPoint (i, pt);
    for (j = 0; j < 3; ++j){
      diff[j] = pt[j] - centre[j];
    }
    // What is the component along the normal?
    compNorm = vtkMath::Dot(diff, normal);

    // Subtract two lots of this component to obtain
    // the reflection.
    for (j = 0; j < 3; ++j){
      pt[j] -= (2 * compNorm * normal[j]);
      // copy.
      ptTemp[j] = pt[j];
    }

    // Now find the nearest point on the surface
    _locator->FindClosestPoint(pt);

    // How far is it from the reflected point?
    sumDist += vtkMath::Distance2BetweenPoints(pt, ptTemp);
  }

//  cout << "bla" << endl;

  return sumDist / (float) noOfPoints;

}


void phiTheta2normal(float phi, float theta, float* normal){
  // First need to convert to radians.
  theta = theta * M_PI / 180.0;
  phi   = phi   * M_PI / 180.0;

  normal[0] = cos(phi) * sin(theta);
  normal[1] = sin(phi) * sin(theta);
  normal[2] = cos(theta);
}

// Work in degrees for easier optimisation
void normal2phiTheta(float *normal, float &phi, float &theta){

  // Add pi to ensure range is in 0->2 pi
  phi = atan2f(normal[1], normal[0]);// + M_PI;

  if (phi < 0)
    phi += 2 * M_PI;

  // r is radius in xy plane
  double r = fabs(sqrt(normal[0]*normal[0] + normal[1]*normal[1]));

  // Adding pi/2 ensures range is from 0 to pi, angle measured from north
  // pole.
  theta = atan2f(r, normal[2]);

  if (theta < 0)
    theta += M_PI_2;

  theta = 180.0 * theta / M_PI;
  phi   = 180.0 * phi / M_PI;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
