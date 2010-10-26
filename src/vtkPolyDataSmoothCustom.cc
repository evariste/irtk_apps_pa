
#include <vtkPolyDataSmoothCustom.h>

vtkPolyDataSmoothCustom::vtkPolyDataSmoothCustom()
{
  _input = NULL;
  _output = NULL;
  _normals = NULL;
  _distances = NULL;

  _CofG[0] = _CofG[1] = _CofG[2] = 0;
}

vtkPolyDataSmoothCustom::~vtkPolyDataSmoothCustom()
{

}

void vtkPolyDataSmoothCustom::Initialize(vtkPolyData *polydata)
{
  _TrackingOn = false;

  _SmoothnessThreshold = -1.0f;

  _NoOfIterations = 1;


  // Normals need to be recalculated before using.
  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SplittingOff();
  normalsFilter->SetInput(polydata);
  normalsFilter->Modified();
  normalsFilter->Update();

  if (_input != NULL)
    _input->Delete();

  if (_output != NULL)
    _output->Delete();

  if (_normals != NULL)
    _normals->Delete();

  if (_distances != NULL)
    _distances->Delete();


  _input = vtkPolyData::New();
  _input = normalsFilter->GetOutput();
  _input->Update();
  _input->BuildCells();
  _input->BuildLinks();

  int noOfPoints = _input->GetNumberOfPoints();

  _normals = vtkFloatArray::New();
  _normals->SetNumberOfComponents(3);
  _normals->SetNumberOfTuples(noOfPoints);
  _normals = (vtkFloatArray*) _input->GetPointData()->GetNormals();


  _distances = vtkFloatArray::New();
  _distances->SetNumberOfComponents(1);
  _distances->SetNumberOfTuples(noOfPoints);
  for (int j = 0; j < noOfPoints; ++j){
    _distances->SetTuple1(j, 0.0);
  }

//  normalsFilter->Delete();
}

void vtkPolyDataSmoothCustom::SetInput(vtkPolyData *polydata)
{
  this->Initialize(polydata);
}

void vtkPolyDataSmoothCustom::Run()
{
  int i, j, k, noOfPoints;
  double E_H2, area, h2norm = 0, totalArea, update[3], currPos[3];
  double v1[3], v2[3], v3[3], centre[3], triangleArea, dx, dy, dz;
  double dist, val, *normal, radiusOld, radiusNew, scaleFactor, cogOld[3], shift[3];
  unsigned short noOfCells = 0;
  vtkIdType* cells = NULL;
  vtkTriangle* triangle = NULL;
  vtkIdList* ptIds = NULL;



  noOfPoints = _input->GetNumberOfPoints();

  // the points array
  double *pts = new double[3*noOfPoints];
  vtkPoints* pts_original = _input->GetPoints();

  for (i = 0; i <= _NoOfIterations; ++i){

    cout << "iteration  " << i << " "; cout.flush();

    // Estimate \int H^2 dA by multiplying E(H^2) with Area.
    E_H2 = this->HSquareRobustMean();
    cout << "bla " << endl;

    area = this->SurfaceArea();


    // The L_2 norm using the Tosun formulation (MedIA 2004)
    h2norm = sqrt(E_H2 * area / 4.0 / M_PI);

    if (h2norm < _SmoothnessThreshold){
      break;
    }

    // Calculate CofG and radius before starting.
    this->UpdateCentreOfGravity();
    cogOld[0] = this->_CofG[0];
    cogOld[1] = this->_CofG[1];
    cogOld[2] = this->_CofG[2];
    radiusOld = this->GetRadius();


    cout << h2norm << endl;

    // Loop over surface.
    for (j = 0; j < noOfPoints; ++j){

      // Initialisation for current point.
      totalArea = 0;
      cells = NULL;

      update[0] = 0;
      update[1] = 0;
      update[2] = 0;

      // Store the current position of the node.
      _input->GetPoint(j, currPos);

      // What cells does this node adjoin?
      _input->GetPointCells(j, noOfCells, cells);

      if ( cells == NULL )
        continue;

      for (k = 0; k < noOfCells; ++k){
        triangle = vtkTriangle::SafeDownCast(_input->GetCell(cells[k]));

        if ( triangle != NULL ){
          ptIds = triangle->GetPointIds();

          _input->GetPoint(ptIds->GetId(0), v1);
          _input->GetPoint(ptIds->GetId(1), v2);
          _input->GetPoint(ptIds->GetId(2), v3);

          triangleArea = vtkTriangle::TriangleArea(v1, v2, v3);
          vtkTriangle::TriangleCenter(v1, v2, v3, centre);

          totalArea += triangleArea;

          update[0] += triangleArea * centre[0];
          update[1] += triangleArea * centre[1];
          update[2] += triangleArea * centre[2];
        }
      }

      update[0] /= totalArea;
      update[1] /= totalArea;
      update[2] /= totalArea;

      dx = _RelaxationFactor * (update[0] - currPos[0]);
      dy = _RelaxationFactor * (update[1] - currPos[1]);
      dz = _RelaxationFactor * (update[2] - currPos[2]);


      pts[j*3]   = currPos[0] + dx;
      pts[j*3+1] = currPos[1] + dy;
      pts[j*3+2] = currPos[2] + dz;

      if (_TrackingOn == true){
        dist = sqrt(dx*dx + dy*dy + dz*dz);

        normal = _normals->GetTuple3(j);
        val = normal[0]*dx + normal[1]*dy + normal[2]*dz;
        if (val < 0){
          dist = -1.0 * dist;
        }

        val = _distances->GetTuple1(j);
        _distances->SetTuple1(j, val + dist);
      }

    }

    for (j = 0; j < noOfPoints; ++j){
      pts_original->SetPoint(j, pts + j*3);
    }

    _input->SetPoints(pts_original);
    _input->Update();



    // update radius and centre of gravity
    this->UpdateCentreOfGravity();
    radiusNew = this->GetRadius();

    shift[0] = cogOld[0] - _CofG[0];
    shift[1] = cogOld[1] - _CofG[1];
    shift[2] = cogOld[2] - _CofG[2];

    scaleFactor = radiusOld / radiusNew;
    this->ShiftAndScaleMesh(shift, scaleFactor);

  }

  cout << "Final iteration : " << i << endl;
  cout << "Final L_2 norm of H^2 (threshold) : " << h2norm << " (" << _SmoothnessThreshold << ")" << endl;

  this->Finalize();
}

void vtkPolyDataSmoothCustom::Finalize()
{

  if (_TrackingOn == true){
    _distances->SetName("smoothingDists");
    _input->GetPointData()->AddArray(_distances);
    _input->Update();
  }

  // Normals need to be recalculated before saving.
  cerr << endl << "Recalculating normals" << endl;
  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SplittingOff();
  normalsFilter->SetInput(_input);
  normalsFilter->Modified();
  normalsFilter->Update();

  _output = vtkPolyData::New();
  _output = normalsFilter->GetOutput();
  _output->Update();

}

vtkPolyData *vtkPolyDataSmoothCustom::GetOutput()
{
  return _output;

}


double vtkPolyDataSmoothCustom::HSquareRobustMean()
{
  int i, noOfPoints;
  double meanValSq = 0.0;
  double val;

  noOfPoints = _input->GetNumberOfPoints();

  vtkCurvatures *curve = vtkCurvatures::New();
  curve->SetInput(_input);
  curve->SetCurvatureTypeToMean();
  curve->Update();

  vtkPolyData *curveOut = vtkPolyData::New();
  curveOut = curve->GetOutput();

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars = (vtkFloatArray *) curveOut->GetPointData()->GetScalars("Mean_Curvature");

  // Get a robust mean.

  double minVal, maxVal;
  double lo = 5.0;
  double hi = 95.0;
  float *data = new float[1 + noOfPoints];
  int count;

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    data[1 + i] = val;
  }

  sort(noOfPoints, data);

  i = 1 + (int) round( (double) lo * (noOfPoints - 1) / 100.0);
  minVal = data[i];
  i = 1 + (int) round( (double) hi * (noOfPoints - 1) / 100.0);
  maxVal = data[i];

  count = 0;
  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    if (val < minVal || val > maxVal)
      continue;

    meanValSq += val * val;
    ++count;
  }

  delete [] data;
//  scalars->Delete();
//  curveOut->Delete();
//  curve->Delete();

  return meanValSq / count;

}


double vtkPolyDataSmoothCustom::SurfaceArea()
{
  vtkCellArray* facets = _input->GetPolys();
  vtkTriangle* facet = vtkTriangle::New();

  double A = 0.0;
  double v0[3], v1[3], v2[3];

  vtkIdType f, *vert = 0;

  facets->InitTraversal();
  while (facets->GetNextCell(f, vert)){

    _input->GetPoint(vert[0], v0);
    _input->GetPoint(vert[1], v1);
    _input->GetPoint(vert[2], v2);

    A += double(facet->TriangleArea(v0, v1, v2));
  }
  facet->Delete();

  return A;
}

void vtkPolyDataSmoothCustom::UpdateCentreOfGravity()
{

  int i, noOfPoints;
  double cofgx, cofgy, cofgz;
  double point[3];

  cofgx = cofgy = cofgz = 0.0;

  noOfPoints = _input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    _input->GetPoint (i, point);
    cofgx += point[0];
    cofgy += point[1];
    cofgz += point[2];
  }

  if (noOfPoints > 0){
    cofgx /= noOfPoints;
    cofgy /= noOfPoints;
    cofgz /= noOfPoints;
  }

  _CofG[0] = cofgx;
  _CofG[1] = cofgy;
  _CofG[2] = cofgz;

}

void vtkPolyDataSmoothCustom::ShiftAndScaleMesh(double *shift, double factor)
{
  int i, j, noOfPoints;
  double vOld[3], vNew[3];

  noOfPoints = _input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; ++i){
    _input->GetPoint(i, vOld);

    for (j = 0; j < 3; ++j){
      vNew[j] = factor * (vOld[j] + shift[j]);
    }

    _input->GetPoints()->SetPoint(i, vNew);
  }
}
double vtkPolyDataSmoothCustom::GetRadius()
{
  int i, noOfPoints;
  double point[3];
  double r, rSum = 0.0;

  noOfPoints = _input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    _input->GetPoint (i, point);

    r = sqrt((point[0]-_CofG[0])*(point[0]-_CofG[0]) +
             (point[1]-_CofG[1])*(point[1]-_CofG[1]) +
             (point[2]-_CofG[2])*(point[2]-_CofG[2]));

    rSum += r;
  }

  return rSum / noOfPoints;

}
