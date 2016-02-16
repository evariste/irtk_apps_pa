////////////////////////////////////////////////////////////////////////
// polydatasmooth with an option to stop when the L2 norm of H (extrinsic
// curvature reaches a certain level indicating a threshold on smoothness..


#if (defined HAS_VTK)

#include <irtkImage.h>
//#include <nr.h>
#include <gsl/gsl_vector.h> /*For Vectors*/
#include <gsl/gsl_sort_vector.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCurvatures.h>

char *input_name = NULL;
char *output_name = NULL;

void usage()
{
  cerr << "" << endl;
  cerr << " Usage: polydatasmooth [input] [output] [iterations] [relaxationFactor]  <options>" << endl;
  cerr << "" << endl;
  cerr << " Area weighted Laplacian relaxation of a surface." << endl;
  cerr << "" << endl;
  cerr << " Iteratively move mesh nodes towards the centroid of the adjacent nodes." << endl;
  cerr << " The relaxation factor (RF) determines how far each node moves towards the " << endl;
  cerr << " local centroid (0<= RF <= 1).  The new position of a node is a weighted " << endl;
  cerr << " combination of its previous position and the neighbours' centroid.  " << endl;
  cerr << " RF = 1 moves a node all the way to the centroid while RF = 0 keeps a " << endl;
  cerr << " node at its previous position." << endl;
  cerr << "" << endl;
  cerr << "Options: " << endl;
  cerr << "" << endl;
  cerr << "-track : " << endl;
  cerr << "         Track the signed distances traversed by points (out is positive)." << endl;
  cerr << "         Store the results in scalar array \"smoothingDists\"" << endl;
  cerr << "-threshold [value] " << endl;
  cerr << "         A value indicating the smoothness as the L_2 norm of H^2 over the surface." << endl;
  cerr << "         See Tosun, MedIA, 2004.  Iterations stop if the norm for the surface drops" << endl;
  cerr << "         drops below the given value." << endl;
  cerr << "" << endl;
  cerr << "" << endl;

  exit(1);
}

double hSquareRobustMean(vtkPolyData* surface){

  int i, noOfPoints;
  double meanVal = 0.0;
  double val;

  noOfPoints = surface->GetNumberOfPoints();

  vtkCurvatures *curve = vtkCurvatures::New();
  curve->SetInputData(surface);
  curve->SetCurvatureTypeToMean();
  curve->Update();

  vtkPolyData *output = vtkPolyData::New();
  output = curve->GetOutput();

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars = (vtkFloatArray *) output->GetPointData()->GetScalars("Mean_Curvature");

  // Get a robust mean.

  double minVal, maxVal;
  double lo = 5.0;
  double hi = 95.0;
  gsl_vector * data = gsl_vector_alloc (noOfPoints);
  //float *data = new float[1 + noOfPoints];
  int count;

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    //data[1 + i] = val;
	gsl_vector_set(data, i, val);
  }

  //sort(noOfPoints, data);
  gsl_sort_vector(data);

  //i = 1 + (int) round( (double) lo * (noOfPoints - 1) / 100.0);
  //minVal = data[i];
  //i = 1 + (int) round( (double) hi * (noOfPoints - 1) / 100.0);
  //maxVal = data[i];

  i = (int) round( (double) lo * (noOfPoints - 1) / 100.0);
  minVal = gsl_vector_get(data, i);
  i = (int) round( (double) hi * (noOfPoints - 1) / 100.0);
  maxVal = gsl_vector_get(data, i);

  count = 0;
  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    if (val < minVal || val > maxVal)
      continue;

    meanVal += val * val;
    ++count;
  }

  return meanVal / count;

}

double surfaceArea(vtkPolyData* surface){
  vtkCellArray* facets = surface->GetPolys();
  vtkTriangle* facet = vtkTriangle::New();

  double A = 0.0;
  double v0[3], v1[3], v2[3];

  vtkIdType f, *vert=0;
  facets->InitTraversal();
  while (facets->GetNextCell(f,vert)){

    surface->GetPoint(vert[0],v0);
    surface->GetPoint(vert[1],v1);
    surface->GetPoint(vert[2],v2);

    A += double(facet->TriangleArea(v0,v1,v2));
  }

  return A;
}

void getCoG(vtkPolyData *input, double*cog){

  int i, noOfPoints;
  double cofgx, cofgy, cofgz;
  double point[3];

  cofgx = cofgy = cofgz = 0.0;

  noOfPoints = input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    input->GetPoint (i, point);
    cofgx += point[0];
    cofgy += point[1];
    cofgz += point[2];
  }

  if (noOfPoints > 0){
    cofgx /= noOfPoints;
    cofgy /= noOfPoints;
    cofgz /= noOfPoints;
  }

  cog[0] = cofgx;
  cog[1] = cofgy;
  cog[2] = cofgz;

}

void shiftAndScalePolyData(vtkPolyData* input, double *shift, double factor){
  int i, j, noOfPoints;
  noOfPoints = input->GetNumberOfPoints();
  double vOld[3], vNew[3];

  for (i = 0; i < noOfPoints; ++i){
    input->GetPoint(i, vOld);

    for (j = 0; j < 3; ++j){
      vNew[j] = factor * (vOld[j] + shift[j]);
    }

    input->GetPoints()->SetPoint(i, vNew);
  }
}

double meanRadius(vtkPolyData* input, double*cog){
  int i, noOfPoints;
  double point[3];
  double r, rSum = 0.0;

  noOfPoints = input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    input->GetPoint (i, point);
    r = sqrt((point[0]-cog[0])*(point[0]-cog[0]) +
             (point[1]-cog[1])*(point[1]-cog[1]) +
             (point[2]-cog[2])*(point[2]-cog[2]));
    rSum += r;
  }

  return rSum / noOfPoints;
}

int main(int argc, char **argv)
{
  int i, j, k;
  bool ok;
  int noOfIterations;
  double relaxationFactor;
  double* pts;
  double E_H2, area;

  double currPos[3];
  unsigned short noOfCells = 0;
  vtkIdType* cells = NULL;
  vtkTriangle* triangle = NULL;
  double totalArea = 0;
  double update[3];
  vtkIdList* ptIds = NULL;
  double v1[3], v2[3], v3[3], centre[3];
  double triangleArea = 0;
  double dx, dy, dz;
  double dist, val;
  double *normal;
  double h2norm;
  double smoothnessThreshold = -1.0f;

  double cogOld[3];
  double cogNew[3];
  double radiusOld, radiusNew;
  double shift[3];
  double scaleFactor;

  int trackingOn = false;

  if (argc < 4){
    usage();
  }

  // Parse arguments.
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;
  noOfIterations = atoi(argv[1]);
  argv++;
  argc--;
  relaxationFactor = atof(argv[1]);
  argv++;
  argc--;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-track") == 0)){
      argc--;
      argv++;
      trackingOn = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;
      smoothnessThreshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
     if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Input        : " << input_name << endl;
  cerr << "Output       : " << output_name << endl;
  cerr << "Iterations   : " << noOfIterations << endl;
  cerr << "Relax factor : " << relaxationFactor << endl;

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();
  vtkPolyData* inputRead = reader->GetOutput();


  // Normals need to be recalculated before using.
  vtkPolyDataNormals *normalsFilter2 = vtkPolyDataNormals::New();
  normalsFilter2->SplittingOff();
  normalsFilter2->SetInputData(inputRead);
  normalsFilter2->Modified();
  normalsFilter2->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = normalsFilter2->GetOutput();


  input->BuildCells();
  input->BuildLinks();

  int noOfPoints = input->GetNumberOfPoints();

  vtkFloatArray *normals = vtkFloatArray::New();
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(noOfPoints);
  normals = (vtkFloatArray*) input->GetPointData()->GetNormals();

  // the points array
  pts = new double[3*noOfPoints];

  vtkFloatArray *dists = vtkFloatArray::New();
  dists->SetNumberOfComponents(1);
  dists->SetNumberOfTuples(noOfPoints);
  for (j = 0; j < noOfPoints; ++j){
    dists->SetTuple1(j, 0.0);
  }

  // perform the smoothing
  vtkPoints* pts_original = input->GetPoints();

  // Iteration loop.

  for (i = 0; i <= noOfIterations; ++i){

    cout << "iteration  " << i << " ";

    // Estimate \int H^2 dA by multiplying E(H^2) with Area.
    E_H2 = hSquareRobustMean(input);
    area = surfaceArea(input);

    // The L_2 norm using the Tosun formulation (MedIA 2004)
    h2norm = sqrt(E_H2 * area / 4.0 / M_PI);

    if (h2norm < smoothnessThreshold){
      break;
    }

    getCoG(input, cogOld);
    radiusOld = meanRadius(input, cogOld);

    cout << h2norm << endl;
//    cout << radiusOld << " / ";
//    cout << cogOld[0] << " " << cogOld[1] << " " << cogOld[2] << endl;

    // Loop over surface.
    for (j = 0; j < noOfPoints; ++j){

      // Initialisation for current point.
      totalArea = 0;
      cells = NULL;

      update[0] = 0;
      update[1] = 0;
      update[2] = 0;

      // Store the current position of the node.
      input->GetPoint(j, currPos);

      // What cells does this node adjoin?
      input->GetPointCells(j, noOfCells, cells);

      if ( cells == NULL )
        continue;

      for (k = 0; k < noOfCells; ++k){
        triangle = vtkTriangle::SafeDownCast(input->GetCell(cells[k]));

        if ( triangle != NULL ){
          ptIds = triangle->GetPointIds();

          input->GetPoint(ptIds->GetId(0), v1);
          input->GetPoint(ptIds->GetId(1), v2);
          input->GetPoint(ptIds->GetId(2), v3);

          triangleArea = vtkTriangle::TriangleArea(v1, v2, v3);
          vtkTriangle::TriangleCenter(v1, v2, v3, centre);

          totalArea += triangleArea;

          update[0] += triangleArea * centre[0];
          update[1] += triangleArea * centre[1];
          update[2] += triangleArea * centre[2];
        }
      }

      if (totalArea <= 0.0){
        update[0] = currPos[0];
        update[1] = currPos[1];
        update[2] = currPos[2];
      } else {
      	update[0] /= totalArea;
      	update[1] /= totalArea;
      	update[2] /= totalArea;
      }

      dx = relaxationFactor * (update[0] - currPos[0]);
      dy = relaxationFactor * (update[1] - currPos[1]);
      dz = relaxationFactor * (update[2] - currPos[2]);


      pts[j*3]   = currPos[0] + dx;
      pts[j*3+1] = currPos[1] + dy;
      pts[j*3+2] = currPos[2] + dz;

      if (trackingOn == true){
        dist = sqrt(dx*dx + dy*dy + dz*dz);

        normal = normals->GetTuple3(j);
        val = normal[0]*dx + normal[1]*dy + normal[2]*dz;
        if (val < 0){
          dist = -1.0 * dist;
        }

        val = dists->GetTuple1(j);
        dists->SetTuple1(j, val + dist);
      }

    }

    for (j = 0; j < noOfPoints; ++j){
      pts_original->SetPoint(j, pts + j*3);
    }

    input->SetPoints(pts_original);


    // update radius and centre of gravity
    getCoG(input, cogNew);
    radiusNew = meanRadius(input, cogNew);

    shift[0] = cogOld[0] - cogNew[0];
    shift[1] = cogOld[1] - cogNew[1];
    shift[2] = cogOld[2] - cogNew[2];

    scaleFactor = radiusOld / radiusNew;

    shiftAndScalePolyData(input, shift, scaleFactor);

//    getCoG(input, cogNew);
//    radiusNew = meanRadius(input, cogNew);
//    cout << " adjusted rad and cog : " <<  radiusNew << " / ";
//    cout << cogNew[0] << " " << cogNew[1] << " " << cogNew[2] << endl;

  }

  cout << "Final iterations : " << i << endl;
  cout << "Final L_2 norm of H^2 (threshold) : " << h2norm << " (" << smoothnessThreshold << ")" << endl;


  if (trackingOn == true){
    dists->SetName("smoothingDists");
    input->GetPointData()->AddArray(dists);
  }

  // Normals need to be recalculated before saving.
  cerr << endl << "Recalculating normals" << endl;
  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SplittingOff();
  normalsFilter->SetInputData(input);
  normalsFilter->Modified();
  normalsFilter->Update();


  // save as a vtk file
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(normalsFilter->GetOutput());
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Write();

  reader->Delete();
  writer->Delete();
  delete [] pts;
  return 0;
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

