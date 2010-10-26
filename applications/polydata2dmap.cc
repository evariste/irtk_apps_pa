#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
// Find the average distance of a set of points from a surface
//

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <irtkLocator.h>

#include <irtkTransformation.h>

char *target_landmarks_name = NULL;
char *source_surface_name = NULL;
char *dofin_name = NULL;

int main(int argc, char **argv){

  bool ok;
  int locatorType;
  irtkTransformation *transformation = NULL;
  bool matching = false;

  if (argc < 3){
    cerr << argv[0] << " [targetLandmarks] [sourceSurface] " << endl;
    cerr << "Find the average distance of a set of points from a surface" << endl;
    cerr << " <-locator>   Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
    cerr << " <-dofin>     transformation file" << endl;
    cerr << " <-matching>  The two surfaces are assumed to have the same connectivity, find mean distance between" << endl;
    cerr << "              corresponding points." << endl;
    exit(1);
  }

  target_landmarks_name = argv[1];
  argv++;
  argc--;
  source_surface_name = argv[1];
  argv++;
  argc--;

  locatorType = 1;

  while (argc > 1){
    ok = false;

    if ((ok == false) && (strcmp(argv[1], "-locator") == 0)){
      argc--;
      argv++;
      locatorType = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-matching") == 0)){
      argc--;
      argv++;
      matching = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  cout << "Reading data ... "; cout.flush();
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(target_landmarks_name);
  reader->Update();

  vtkPolyData* targetLandmarks = vtkPolyData::New();
  targetLandmarks = reader->GetOutput();

  vtkPolyDataReader *reader2 = vtkPolyDataReader::New();
  reader2->SetFileName(source_surface_name);
  reader2->Update();

  vtkPolyData* sourceSurface = vtkPolyData::New();
  sourceSurface = reader2->GetOutput();
  cout << "done" << endl;

  // Create locator
  irtkLocator *source_locator = new irtkLocator;
  source_locator->SelectLocatorType(locatorType);
  source_locator->SetDataSet(sourceSurface);

  // Create transformation
  if (dofin_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dofin_name);
  } else {
    transformation = new irtkRigidTransformation;
  }

  int n = targetLandmarks->GetNumberOfPoints();

  double target_point[3], tmp_point[3], source_point[3];
  double error = 0.0;
  double squaredDist;
  double sumSquaredDists = 0.0;

  if (matching == false){
    for (int i = 0; i < n; i++){
      targetLandmarks->GetPoints()->GetPoint (i, target_point);
      tmp_point[0] = target_point[0];
      tmp_point[1] = target_point[1];
      tmp_point[2] = target_point[2];
      transformation->Transform(tmp_point[0], tmp_point[1], tmp_point[2]);
      source_point[0] = tmp_point[0];
      source_point[1] = tmp_point[1];
      source_point[2] = tmp_point[2];
      (void) source_locator->FindClosestPoint (source_point);

      squaredDist = (tmp_point[0] - source_point[0]) * (tmp_point[0] - source_point[0]) +
                    (tmp_point[1] - source_point[1]) * (tmp_point[1] - source_point[1]) +
                    (tmp_point[2] - source_point[2]) * (tmp_point[2] - source_point[2]);
      if (fabs(squaredDist) > FLT_MIN){
        error += sqrt(squaredDist);
      }
    }
  } else {

    if (n != sourceSurface->GetNumberOfPoints()){
      cerr << "Polydata sets must have same number of points if -match option selected." << endl;
      exit(1);
    }

    for (int i = 0; i < n; i++){
      targetLandmarks->GetPoints()->GetPoint (i, target_point);
      transformation->Transform(target_point[0], target_point[1], target_point[2]);

      sourceSurface->GetPoints()->GetPoint (i, source_point);

      squaredDist = (target_point[0] - source_point[0]) * (target_point[0] - source_point[0]) +
                    (target_point[1] - source_point[1]) * (target_point[1] - source_point[1]) +
                    (target_point[2] - source_point[2]) * (target_point[2] - source_point[2]);

      if (fabs(squaredDist) > FLT_MIN){
        error += sqrt(squaredDist);
        sumSquaredDists += squaredDist;
      }
    }

    cout << "Mean Squared Dist: " << sumSquaredDists / ((double) n) << endl;
    cout << "R.M.S.           : " << sqrt(sumSquaredDists / ((double) n)) << endl;
  }

  cout << "Number of points : " << n << endl;
  cout << "Average distance : " << error / ((double) n) << endl;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
