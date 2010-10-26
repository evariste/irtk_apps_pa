#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
// Evolve a surface outwards to a cortical boundary in a T1 image
//

#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>

#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>

#include <vtkPointLocator.h>
#include <irtkImage.h>
#include <irtkImageFunction.h>

char *input_surface_name = NULL;
char *input_image_name   = NULL;
char *output_image_name  = NULL;
 
void WriteSurfaceToImage(vtkPolyData *surf, irtkGreyImage *image);

int expandSurface(vtkPolyData*,vtkDataArray*, vtkPointLocator*, irtkGreyImage*, irtkInterpolateImageFunction*, double, int, int, double, double);

void fillCenter(irtkGreyImage *);
void Label3D(int, int, int, irtkGreyImage *);

void GetStatistics(vtkPolyData *, irtkGreyImage *, irtkInterpolateImageFunction *, double &, double &);

void usage(char **argv){
  cout << argv[0] << " [imageIn] [surfaceIn] [imageOut] [stepLength]" << endl;
  exit(1);
}

int main(int argc, char **argv){

  int i, iterations;
  int noOfSteps = 50;
  double stepLength = 0.1;
  int neighbourhood = 10;
  double mean, stdev;
  iterations = 6;

  irtkGreyImage *input         = NULL;
  vtkPolyData* surface        = NULL;
  vtkDataArray *normals       = NULL;

  if (argc < 4){
    usage(argv);
  }

  //
  input_image_name = argv[1];
  argv++;
  argc--;
  input_surface_name = argv[1];
  argv++;
  argc--;
  output_image_name = argv[1];
  argv++;
  argc--;
  stepLength = atof(argv[1]);
  argv++;
  argc--;

  //
  input = new irtkGreyImage;
  input->Read(input_image_name);

  irtkInterpolateImageFunction *interpolator = NULL;
  interpolator = new irtkLinearInterpolateImageFunction;
  interpolator->SetInput(input);
  interpolator->Initialize();

  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_surface_name);
  reader->Update();

  surface = vtkPolyData::New();
  surface = reader->GetOutput();

  normals = vtkFloatArray::New();
  normals->SetNumberOfComponents(3);
  normals = surface->GetPointData()->GetNormals();

  vtkPointLocator *locator = vtkPointLocator::New();
  locator->SetDataSet(surface);
  locator->BuildLocator();

  vtkPolyDataConnectivityFilter *connFilter = vtkPolyDataConnectivityFilter::New();
  connFilter->SetExtractionModeToLargestRegion();
  connFilter->SetInput(surface);
  connFilter->Modified();
  connFilter->Update();
  surface = connFilter->GetOutput();
  surface->Modified();
  surface->Update();

  vtkSmoothPolyDataFilter *smoother = vtkSmoothPolyDataFilter::New();
  smoother->SetNumberOfIterations(100);
  smoother->BoundarySmoothingOn();
  smoother->SetFeatureAngle(120);
  smoother->SetEdgeAngle(90);
  smoother->SetRelaxationFactor(.025);

  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SplittingOff();

  //
  GetStatistics(surface, input, interpolator, mean, stdev);

  //


  // noOfSteps now redundant.

  for (i = 0; i < iterations; ++i){

    expandSurface(surface, normals, locator, input, interpolator, 
                  stepLength, neighbourhood, noOfSteps, mean, stdev);

    stepLength *= 0.7;

  }


    smoother->SetInput(surface);
    smoother->Modified();
    smoother->Update();
    surface = smoother->GetOutput();
    surface->Modified();
    surface->Update();

    normalsFilter->SetInput(surface);
    normalsFilter->Modified();
    normalsFilter->Update();
    surface = normalsFilter->GetOutput();
    surface->Modified();
    surface->Update();


  //  Write surface boundary voxels to an image.
  WriteSurfaceToImage(surface, input);

    fillCenter(input);

  input->Write(output_image_name);

//   vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//   writer->SetInput(surface);
//   writer->SetFileName("modifiedSurf.vtk");
//   writer->Write();


}

void GetStatistics(vtkPolyData *surface, irtkGreyImage *input, irtkInterpolateImageFunction *interpolator, double &mean, double &stdev)
{

  int n, nPts;
  double int_x1, int_y1, int_z1, int_x2, int_y2, int_z2;
  double pt[3];
  double sum, sumsq;
  int count;
  double val;

  nPts = surface->GetNumberOfPoints();
  interpolator->Inside(int_x1, int_y1, int_z1, int_x2, int_y2, int_z2);

  sum = sumsq = 0.0;
  count = 0;


  for (n = 0; n < nPts; ++n){

    (surface->GetPoints())->GetPoint(n, pt);

    input->WorldToImage(pt[0], pt[1] , pt[2]);

    if (pt[0] < int_x1 || pt[0] > int_x2 || pt[1] < int_y1 || pt[1] > int_y2 || pt[2] < int_z1 || pt[2] > int_z2){
      break;
    }

    val = interpolator->EvaluateInside(pt[0], pt[1] , pt[2]);

    sum += val;
    sumsq += val * val;
    ++count;

  }


  mean = sum / ((double) count);
  stdev = sqrt( ( sumsq / ((double) count) ) - mean * mean );

  cout << "GetStatistics :" << endl;
  cout << "count " << count << " nPts " << nPts << endl;
  cout << "mean  " << mean << endl;
  cout << "stdev " << stdev << endl;

}

void fillCenter(irtkGreyImage *input)
{
  // Image should be a binary representation of cortical surface.

  // Find the centre.

  int i, j, k;
  int xdim, ydim, zdim;
  double xmu, ymu, zmu;
  xmu = ymu = zmu = 0;
  int count = 0;

  xdim = input->GetX();
  ydim = input->GetY();
  zdim = input->GetZ();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        if (input->Get(i, j, k) > 0){
          xmu += i;
          ymu += j;
          zmu += k;
          ++count;
        }
      }
    }
  }

  if (count > 0){
    xmu /= count;
    ymu /= count;
    zmu /= count;
  } else {
    cerr << "fillCenter : No voxels labeled > 0 in image." << endl;
    exit(1);
  }

  cout << "CofG at : " << xmu << " " << ymu << " " << zmu << endl;
  i = round(xmu);
  j = round(ymu);
  k = round(zmu);

  Label3D(i, j, k, input);


}


void Label3D(int x, int y, int z, irtkGreyImage *input)
{
  if ((x >= 0) && (x < input->GetX()) && (y >= 0) && (y < input->GetY()) && 
      (z >= 0) && (z < input->GetZ())){
    if (input->Get(x, y, z) != 1){
      input->Put(x, y, z, 1);

      Label3D(x-1, y  , z  , input);
      Label3D(x+1, y  , z  , input);
      Label3D(x  , y-1, z  , input);
      Label3D(x  , y+1, z  , input);
      Label3D(x  , y  , z-1, input);
      Label3D(x  , y  , z+1, input);
    }
  }
}



int expandSurface(vtkPolyData *surface, vtkDataArray *normals, vtkPointLocator *locator, irtkGreyImage *input,
                  irtkInterpolateImageFunction *interpolator,
                  double stepLength, int neighbourhood, int noOfSteps, double mean, double stdev)
{
  int i, n, nPts;
  bool ok;
  int nearestPtCount;
  int numberModified = 0;
  double pt[3];
  double nml[3], nmlSum[3];
  double mag;
  double u, x, y, z, xold, yold, zold, xnew, ynew, znew, dx, dy, dz;
  double int_x1, int_y1, int_z1, int_x2, int_y2, int_z2;
  double intensity;
  double *profile = new double[noOfSteps];
  double upperLimit, lowerLimit;
  double baseline;

  upperLimit = mean + 2.5 * stdev;
  lowerLimit = mean - 2.5 * stdev;


  vtkIdList *nearestPtIDs = vtkIdList::New();

  nPts = surface->GetNumberOfPoints();

  vtkDataArray *displacements;
  displacements = vtkFloatArray::New();
  displacements->SetNumberOfComponents(3);
  displacements->SetNumberOfTuples(nPts);

  interpolator->Inside(int_x1, int_y1, int_z1, int_x2, int_y2, int_z2);

  for (n = 0; n < nPts; ++n){

    (surface->GetPoints())->GetPoint(n, pt);

    locator->FindClosestNPoints(neighbourhood, pt, nearestPtIDs);
    nearestPtCount = nearestPtIDs->GetNumberOfIds();

    nmlSum[0] = nmlSum[1] = nmlSum[2] = 0;    

    if (nearestPtCount > 0){
      for (i = 0; i < nearestPtCount; ++i){
        normals->GetTuple(nearestPtIDs->GetId(i), nml);
        nmlSum[0] += nml[0];
        nmlSum[1] += nml[1];
        nmlSum[2] += nml[2];
      }
      nml[0] = nmlSum[0] / nearestPtCount;
      nml[1] = nmlSum[1] / nearestPtCount;
      nml[2] = nmlSum[2] / nearestPtCount;
    } else {
      normals->GetTuple(n, nml);
    }

    mag = sqrt(nml[0] * nml[0] + nml[1] * nml[1] + nml[2] * nml[2]);
    nml[0] /= mag;
    nml[1] /= mag;
    nml[2] /= mag;


    xold = pt[0];
    yold = pt[1];
    zold = pt[2];
    input->WorldToImage(xold, yold, zold);

    xnew = pt[0] + stepLength * nml[0];
    ynew = pt[1] + stepLength * nml[1];
    znew = pt[2] + stepLength * nml[2];
    input->WorldToImage(xnew, ynew, znew);

    dx = xnew - xold;
    dy = ynew - yold;
    dz = znew - zold;

    x = xold;
    y = yold;
    z = zold;

    if (x < int_x1 || x > int_x2 || y < int_y1 || y > int_y2 || z < int_z1 || z > int_z2){
      continue;
    }

    baseline = interpolator->EvaluateInside(x, y, z);

    if (baseline < lowerLimit || baseline > upperLimit){
      continue;
    }

    ok = true;

    for (u = 0.0; u < 1.0; u += 0.1){
      x = xold + u * dx;
      y = yold + u * dy;
      z = zold + u * dz;

      if (x < int_x1 || x > int_x2 || y < int_y1 || y > int_y2 || z < int_z1 || z > int_z2){
        ok = false;
        break;
      }

      intensity = interpolator->EvaluateInside(x, y, z);

      if (fabs(intensity - baseline) > 0.2 * baseline){
        ok = false;
        break;
      }
//       if (intensity < lowerLimit || intensity > upperLimit){
//         ok = false;
//         break;
//       }
    }

    if (ok == false)
      continue;


//     if (xold < int_x1 || xold > int_x2 || yold < int_y1 || yold > int_y2 || zold < int_z1 || zold > int_z2){
//       continue;
//     }

//     intensity = interpolator->EvaluateInside(xold, yold, zold);

//     if (intensity < lowerLimit || intensity > upperLimit)
//       continue;

//     if (xnew < int_x1 || xnew > int_x2 || ynew < int_y1 || ynew > int_y2 || znew < int_z1 || znew > int_z2){
//       continue;
//     }

//     intensity = interpolator->EvaluateInside(xnew, ynew, znew);

//     if (intensity < lowerLimit || intensity > upperLimit)
//       continue;




    pt[0] += stepLength * nml[0];
    pt[1] += stepLength * nml[1];
    pt[2] += stepLength * nml[2];

    (surface->GetPoints())->SetPoint(n, pt);

    ++numberModified;


//     dx = x2 - x1;
//     dy = y2 - y1;
//     dz = z2 - z1;

//     for (step = 0; step < noOfSteps; ++step){
//       profile[step] = -1;
//     }

//     intensity = interpolator->EvaluateInside(x1, y1, z1);

//     for (step = 0; step < noOfSteps; ++step){
//       x = x1 + dx * step;
//       y = y1 + dy * step;
//       z = z1 + dz * step;

//       if (x < int_x1 || x > int_x2 || y < int_y1 || y > int_y2 || z < int_z1 || z > int_z2){
//         cout << "xx " << endl;
//         break;
//       }

//       profile[step] = interpolator->EvaluateInside(x, y, z); 

//     }

//     step = 0;
//     while (step < noOfSteps && profile[step] / intensity > 0.8){
//       step++;
//     }

//     if (step >= noOfSteps)
//       step = 0;

//     pt[0] = x1 + dx * step;
//     pt[1] = y1 + dy * step;
//     pt[2] = z1 + dz * step;
//     input->ImageToWorld(pt[0], pt[1], pt[2]);
//     (surface->GetPoints())->SetPoint(n, pt);

  }

  surface->Modified();
  surface->Update();


  delete [] profile;

  cout << "numberModified  " << numberModified << endl;

  return numberModified;

}


void WriteSurfaceToImage(vtkPolyData *surf, irtkGreyImage *image){

  irtkGreyPixel *ptrPix;
  int voxels, i, j, k;
  int nPoints, nCells, n;
  double x, y, z;
  int xdim, ydim, zdim;
  double a[3], b[3], c[3];
  irtkVector ab, ac, cross;
  double u, v;

  ab.Initialize(3);
  ac.Initialize(3);
  cross.Initialize(3);

  vtkIdList *ptIds = vtkIdList::New();

  ptrPix  = image->GetPointerToVoxels();
  voxels  = image->GetNumberOfVoxels();
  nPoints = surf->GetNumberOfPoints();
  nCells = surf->GetNumberOfCells();

  cout << "nPoints : " << nPoints << endl;

  for (i = 0; i < voxels; ++i){
    *ptrPix = 0;
    ++ptrPix;
  }

  xdim = image->GetX();
  ydim = image->GetY();
  zdim = image->GetZ();

  vtkDataArray *nmls;
  nmls = vtkFloatArray::New();
  nmls->SetNumberOfComponents(3);
  nmls = surf->GetPointData()->GetNormals();

  for (n = 0; n < nCells; ++n){

    surf->GetCellPoints(n, ptIds);

    if (ptIds->GetNumberOfIds() != 3){
      cerr << "More than 3 points in a cell " << endl;
      exit(1);
    }

    (surf->GetPoints())->GetPoint(ptIds->GetId(0), a);
    (surf->GetPoints())->GetPoint(ptIds->GetId(1), b);
    (surf->GetPoints())->GetPoint(ptIds->GetId(2), c);

    for (i = 0; i < 3; ++i){
      ab(i) = b[i] - a[i];
      ac(i) = c[i] - a[i];
    }
    cross = ab.CrossProduct(ac);

    // Cover the surface of the current polygon.
    for (u = 0; u < 1; u += 0.1){
      for (v = 0; v < 1; v += 0.1){
//         for (w = 0; w > -1.0; w -= 0.1){

          x = a[0] + (u * ab(0) ) + (v * ac(0) );
          y = a[1] + (u * ab(1) ) + (v * ac(1) );
          z = a[2] + (u * ab(2) ) + (v * ac(2) );

          image->WorldToImage(x, y, z);

          i = (int) round(x);
          j = (int) round(y);
          k = (int) round(z);

          if (i < 0 || i >= xdim || j < 0 || j >= ydim || k < 0 || k >= zdim){
            cerr << "xx " << i << " " << j << " " << k << endl;
            continue;
          }

          image->Put(i, j, k, 1);

//         }
      }
    }
  }

}

//     (surf->GetPoints())->GetPoint(n, pt);

//     nmls->GetTuple(n, nml);

//     x1 = pt[0];
//     y1 = pt[1];
//     z1 = pt[2];

//     image->WorldToImage(x1, y1, z1);

//     x2 = pt[0] - 0.1 * nml[0];
//     y2 = pt[1] - 0.1 * nml[1];
//     z2 = pt[2] - 0.1 * nml[2];

//     image->WorldToImage(x2, y2, z2);

//     dx = x2 - x1;
//     dy = y2 - y1;
//     dz = z2 - z1;

//     for (step = 0; step < 10; ++step){

//       x = x1 + step * dx;
//       y = y1 + step * dy;
//       z = z1 + step * dz;

//       i = (int) round(x);
//       j = (int) round(y);
//       k = (int) round(z);

//       if (i < 0 || i >= xdim || j < 0 || j >= ydim || k < 0 || k >= zdim){
//         cerr << "xx " << i << " " << j << " " << k << endl;
//         continue;
//       }

//       image->Put(i, j, k, 1);


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
