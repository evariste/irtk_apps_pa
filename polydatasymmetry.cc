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
float directionalBounds(float *normal, float *centre);
float objectiveFuncForNR(float params[]);
void writePlaneAsPolydata(float *normal, float *centre);
float distanceToPlane(double x, double y, double z, float *normal, float *centre);

char *surface_name = NULL;
char *template_image_name = NULL;
char *output_image_name = NULL;
char *output_surface_name = NULL;

vtkPolyData* _surface;
irtkLocator *_locator;

void usage()
{
  cerr << " polydatasymmetry [inputSurface] [inputImage] [outputImage] <options>" << endl;
  cerr << " " << endl;
  cerr << " Estimate the plane of symmetry of [inputSurface]." << endl;
  cerr << " Write result as a signed distance map to [outputImage]." << endl;
  cerr << " [inputImage] provides the lattice to be used for the output." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -polydata [name.vtk] : Write a visualisation of the plane as a " << endl;
  cerr << "                        vtkPolyData object to the named file." << endl;
  cerr << " -point x y z         : Coordinates of a point that is required to be on" << endl;
  cerr << "                        on the 'positive' side of the plane, i.e. in the" << endl;
  cerr << "                        same direction as the normal.  The calculated normal" << endl;
  cerr << "                        will be flipped to meet this requirement if necessary." << endl;
  cerr << "" << endl;
  exit(0);
}

int main(int argc, char **argv){

  int ok;
  int noOfPoints;
  double cx, cy, cz;
  double pt[3];
  float normal[3], centre[3];
  double xx, xy, xz, yy, yz, zz;
  int i, j, k;
  double x, y, z;
  double val;
  float phi, theta;
  int hasPositivePoint = False;
  float posiPoint[3];
  int rankMeasure, minRankMeasure;
  float minSymmMeasure;
  int minCol = 0;
  float params[6];
  float symmRank[3] = {0, 1, 2};
  float boundsRank[3] = {0, 1, 2};
  float symmMeasures[3];
  float boundsMeasures[3];
  int normalInit = False;



  if (argc < 2){
    usage();
  }

  surface_name = argv[1];
  argv++;
  argc--;
  template_image_name = argv[1];
  argv++;
  argc--;
  output_image_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-polydata") == 0)){
       argc--;
       argv++;
       output_surface_name = argv[1];
       argc--;
       argv++;
       ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-point") == 0)){
        argc--;
        argv++;
        hasPositivePoint = True;
        posiPoint[0] = atof(argv[1]);
        argc--;
        argv++;
        posiPoint[1] = atof(argv[1]);
        argc--;
        argv++;
        posiPoint[2] = atof(argv[1]);
        argc--;
        argv++;
        ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-normal") == 0)){
      argc--;
      argv++;
      normalInit = True;
      normal[0] = atof(argv[1]);
      argc--;
      argv++;
      normal[1] = atof(argv[1]);
      argc--;
      argv++;
      normal[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }

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

  cout << "Centre of gravity (initial estimate for plane centre): " << endl;
  cout << "           (" << cx << ", " << cy << ", " << cz << ")" <<  endl;

  centre[0] = cx;
  centre[1] = cy;
  centre[2] = cz;


  // Set up locator.
  _locator = new irtkLocator;
  _locator->SelectLocatorType(1);
  _locator->SetDataSet(_surface);

  if (normalInit == True){
    cout << "Using initial normal estimate supplied :" << endl;
    cout << "    " << normal[0] << " " << normal[1] << " " << normal[2] << endl;
    val = sqrt(vtkMath::Dot(normal, normal));
    for (i = 0; i < 3; ++i){
      normal[i] = normal[i] / val;
    }
    cout << "Normalised : " << endl;
    cout << "    " << normal[0] << " " << normal[1] << " " << normal[2] << endl;

  } else {


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

    //  cout << "Finding covariance" << endl;

    cov(0, 0) = xx;
    cov(1, 1) = yy;
    cov(2, 2) = zz;

    cov(0, 1) = xy;
    cov(1, 0) = xy;

    cov(0, 2) = xz;
    cov(2, 0) = xz;

    cov(1, 2) = yz;
    cov(2, 1) = yz;

    //  cout << "Finding eigenstuff " << endl;

    irtkMatrix evecs;
    irtkVector evals;


    cov.Eigenvalues(evecs, evals);

    cout << "Covariance matrix for all points:" << endl;
    cov.Print();
    cout << " " << endl;

    cout << "Eigenvectors: " << endl;
    evecs.Print();
    cout << " " << endl;

    cout << "Eigenvalues: " << endl;
    evals.Print();
    cout << " " << endl;

    // Which eigenvector is the best as a plane normal for a symmetry plane?
    minRankMeasure = INT_MAX;
    minSymmMeasure = FLT_MAX;

    cout << "Symmetry and bounds measures and ranks for eigenvectors: " << endl;
    for (j = 0; j < 3; ++j){
      for (i = 0; i < 3; ++i){
        normal[i] = evecs(i, j);
      }

      symmMeasures[j] = symmetryMeasure(normal, centre);
      boundsMeasures[j] = directionalBounds(normal, centre);
    }

    float a[4], b[4];

    for (i = 1; i <= 3; ++i){
      a[i] = symmMeasures[i-1];
      b[i] = i - 1;
    }

    sort2(3, a, b);
    for (i = 1; i <= 3; ++i){
      //    printf("%d\t %2.2f \t %f \n", i, a[i], b[i]);
      symmRank[(int)(round(b[i]))] = i;
    }
    cout << endl;

    for (i = 1; i <= 3; ++i){
      a[i] = boundsMeasures[i-1];
      b[i] = i - 1;
    }
    sort2(3, a, b);
    for (i = 1; i <= 3; ++i){
      //    printf("%d\t %2.2f \t %f \n", i, a[i], b[i]);
      boundsRank[(int)(round(b[i]))] = i;
    }

    cout << endl;
    cout << "----------------" << endl;
    cout << endl;


    cout << "Symmetry:" << endl;
    printf("n\t meas. \t rank\n");
    for (j = 0; j < 3; ++j){
      printf("%d\t %2.2f \t %d \n", j+1, symmMeasures[j], (int)symmRank[j]);
    }
    cout << "Bounds:" << endl;
    printf("n\t meas. \t rank\n");
    for (j = 0; j < 3; ++j){
      printf("%d\t %2.2f \t %d \n", j+1, boundsMeasures[j], (int)boundsRank[j]);
    }
    cout << endl;

    for (j = 0; j < 3; ++j){
      rankMeasure = (int)(symmRank[j] + boundsRank[j]);
      if (rankMeasure == minRankMeasure && symmMeasures[j] < minSymmMeasure){
        minRankMeasure = rankMeasure;
        minSymmMeasure = symmMeasures[j];
        minCol = j;

      }
      if (rankMeasure < minRankMeasure){
        minRankMeasure = rankMeasure;
        minSymmMeasure = symmMeasures[j];
        minCol = j;
      }
    }


    // Assign minimising e-vec to the normal.
    for (i = 0; i < 3; ++i){
      normal[i] = evecs(i, minCol);
    }

    cout << "Eigenvector " << minCol + 1 << " gives minimum combined measure." << endl;
    cout << "Initial normal estimate : " << normal[0] << " " << normal[1] << " " << normal[2] << endl;


  } // else if normalInit == True





  // Express normal in degrees:
  normal2phiTheta(normal, phi, theta);

  // Set up parameters for NR call.
  params[1] = phi;
  params[2] = theta;
  params[3] = centre[0];
  params[4] = centre[1];
  params[5] = centre[2];

  // Bits and pieces needed for call to Powell method.

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

  // Tolerance for stopping condition.  Not really sure what this means.
  float ftol;
  ftol = 0.5;

  // Number of iterations taken.
  int iters;

  // Final value of objective function.
  float fReturned;

  cout << "Powell minimisation ..";
  cout.flush();
  powell(params, initDirs, 5, ftol, &iters, &fReturned, objectiveFuncForNR);
  cout << ". done." << endl;

  phi = params[1];
  theta = params[2];
  centre[0] = params[3];
  centre[1] = params[4];
  centre[2] = params[5];
  phiTheta2normal(phi, theta, normal);

  cout << "Estimates after optimisation: " << endl;
  cout << "  Normal       : " << normal[0] << " " << normal[1] << " " << normal[2] << endl;
  cout << "  Plane centre : " << centre[0] << " " << centre[1] << " " << centre[2] << endl;
  cout << "  Measure      : " << fReturned << endl << endl;


  // Ensure unit normal.
//  val = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  val = sqrt(vtkMath::Dot(normal, normal));
  for (i = 0; i < 3; ++i){
    normal[i] = normal[i] / val;
  }

  // Check if we need to flip the normal.
  if (hasPositivePoint == True){
    for (i = 0; i < 3; ++i){
      posiPoint[i] -= centre[i];
    }
    if (vtkMath::Dot(normal, posiPoint) < 0){
      for (i = 0; i < 3; ++i){
        normal[i] *= -1.0;
      }
      cout << "Normal flipped: " << endl;
      cout << "Normal       : " << normal[0] << " " << normal[1] << " " << normal[2] << endl;
    }
  }

  // Write distance map.
  irtkRealImage *image = new irtkRealImage(template_image_name);

  for (k = 0; k < image->GetZ(); ++k){
    for (j = 0; j < image->GetY(); ++j){
      for (i = 0; i < image->GetX(); ++i){
        x = i;
        y = j;
        z = k;
        image->ImageToWorld(x, y, z);
        val = distanceToPlane(x, y, z, normal, centre);
        image->Put(i, j, k, val);
      }
    }
  }

  image->Write(output_image_name);

  if (output_surface_name != NULL){
    writePlaneAsPolydata(normal, centre);
  }


  for (i = 0; i < 6; ++i){
      delete [] initDirs[i];
  }
  delete [] initDirs;

}

float distanceToPlane(double x, double y, double z, float *normal, float *centre)
{
  x -= centre[0];
  y -= centre[1];
  z -= centre[2];

  return x*normal[0] + y*normal[1] + z*normal[2];
}

// Function assumes that global variable 'output_surface_name'
// is non-null.
void writePlaneAsPolydata(float *normal, float *centre)
{
  // Now make a plane.
  double surfaceBounds[6];
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double nDotC, val;
  int maxInd;
  int i;
  double maxComp;

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

  // Which normal component is largest?
  maxInd = 0;
  maxComp = fabs(normal[0]);
  for (i = 1; i < 3; ++i){
    if (fabs(normal[i]) > maxComp){
      maxInd = i;
      maxComp = fabs(normal[i]);
    }
  }

  nDotC = vtkMath::Dot(normal, centre);

  // Try and get points for the plane that cover the object
  if (maxInd == 2){
    // Normal is mostly in the z direction.
    val = (nDotC - normal[0]*xmin - normal[1]*ymax) / normal[2];
    plane->SetOrigin(xmin, ymax, val);
    val = (nDotC - normal[0]*xmin - normal[1]*ymin) / normal[2];
    plane->SetPoint1(xmin, ymin, val);
    val = (nDotC - normal[0]*xmax - normal[1]*ymax) / normal[2];
    plane->SetPoint2(xmax, ymax, val);
  } else if (maxInd == 1){
    // Normal is mostly in the y direction.
    val = (nDotC - normal[0]*xmin - normal[2]*zmax) / normal[1];
    plane->SetOrigin(xmin, val, zmax);
    val = (nDotC - normal[0]*xmin - normal[2]*zmin) / normal[1];
    plane->SetPoint1(xmin, val, zmin);
    val = (nDotC - normal[0]*xmax - normal[2]*zmax) / normal[1];
    plane->SetPoint2(xmax, val, zmax);
  } else {
    // Normal is mostly in the x direction.
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
  writer->SetFileName(output_surface_name);
  writer->Update();
  writer->Write();

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

//  cout << "  " << phi << " " << theta << " " << centre[0] << " " << centre[1] << " " << centre [2] << " " << measure << endl;

  measure = symmetryMeasure(phi, theta, centre);


  return measure;
}

float symmetryMeasure(float phi, float theta, float *centre)
{
  float normal[3];
  phiTheta2normal(phi, theta, normal);
  return symmetryMeasure(normal, centre);
}

// For a given direction and centre, find the largest extent of the surface along the direction
float directionalBounds(float *normal, float *centre)
{
  int i, j, noOfPoints;
  double pt[3];
  float diff[3];
  double temp, compNorm;
  float minComp, maxComp;

  // Ensure unit normal.
  temp = vtkMath::Norm(normal);
  for (j = 0; j < 3; ++j)
    normal[j] /= temp;

  noOfPoints = _surface->GetNumberOfPoints();

  minComp = FLT_MAX;
  maxComp = -1.0 * FLT_MAX;

  for (i = 0; i < noOfPoints; ++i){
    _surface->GetPoint (i, pt);
    for (j = 0; j < 3; ++j){
      diff[j] = pt[j] - centre[j];
    }

    // What is the component along the normal?
    compNorm = vtkMath::Dot(diff, normal);

    if (maxComp < compNorm)
      maxComp = compNorm;

    if (minComp > compNorm)
      minComp = compNorm;

  }

//   cout << " Min / max comp: " << minComp << " " << maxComp << endl;
//  cout << " Difference    : " << maxComp - minComp << endl;

  return maxComp - minComp;
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
