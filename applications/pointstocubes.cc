#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
// Fatten some points for viewing.

#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCubeSource.h>
#include <vtkSphereSource.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkGlyph3D.h>
#include <irtkImage.h>


char *input_landmarks_name = NULL;
char *output_name = NULL;

int faces[6][4] = {{0, 4, 6, 2}, {1, 3, 7, 5}, {4, 5, 7, 6}, {0, 2, 3, 1}, {2, 6, 7, 3}, {0, 1, 5, 4}};

void usage(){

  cout << "pointstocubes [inputPoints.vtk] [outputPolys.vtk] <-r> <-scalar>" << endl;
  cout << "" << endl;
  cout << "Place a cube centered at each of the input points." << endl;
  cout << "default radius is 1, can assign a scalar value if wanted." << endl;

  exit(1);
}

void putCube(vtkPoints* points, vtkCellArray *polys, int count, double *p, double r){

  int modifiedFaces[6][4];
  double x[8][3];
  int i, j, k, n;

  n = 0;
  for (i = -1; i < 2; i += 2){
    for (j = -1; j < 2; j += 2){
      for (k = -1; k < 2; k += 2){
        x[n][0] = p[0] + i * r;
        x[n][1] = p[1] + j * r;
        x[n][2] = p[2] + k * r;
        n++;
      }
    }
  }

  for (i = 0; i < 8; i++){
    points->InsertPoint(count + i, x[i]);
  }

  for (i = 0; i < 6; i++){

    for (j = 0; j < 4; j++){
      modifiedFaces[i][j] = faces[i][j] + count;
    }

    polys->InsertNextCell(4, modifiedFaces[i]);

  }

} 

int main(int argc, char **argv){

  bool ok;

  if (argc < 3)
    usage();


  input_landmarks_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  double radius = 1.0;
  double scalar = 1.0;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-r") == 0)){
      argc--;
      argv++;
      radius = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-scalar") == 0)){
      argc--;
      argv++;
      scalar = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }



  cout << "Reading target landmarks ... "; cout.flush();
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_landmarks_name);
  reader->Update();
  vtkPolyData* targetLandmarks = vtkPolyData::New();
  targetLandmarks = reader->GetOutput();
  cout << "done" << endl;

  int n = targetLandmarks->GetNumberOfPoints();

  vtkPoints     *points  = vtkPoints::New();
  vtkCellArray  *polys   = vtkCellArray::New();
  vtkFloatArray *scalars = vtkFloatArray::New();

  vtkPolyData   *output  = vtkPolyData::New();
  output->Allocate();

  double p[3];
  int count = 0;

  for (int i = 0; i < n; i++){
    targetLandmarks->GetPoint(i, p);
    // Place a cube centered at p with radius 'radius'.
    putCube(points, polys, count,  &p[0], radius);
    count += 8;
  }

  for (int i = 0; i < count; i++){
    scalars->InsertTuple1(i, scalar);
  }


  output->SetPoints(points);
  output->SetPolys(polys);
  output->GetPointData()->AddArray(scalars);

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(output);
  writer->SetFileName(output_name);
  writer->Write();




  vtkGlyph3D *glyphs = vtkGlyph3D::New();
  
  vtkCubeSource *cube = vtkCubeSource::New();
  cube->SetBounds(-1, 1, -1, 1, -1, 1);

  vtkSphereSource *sphere = vtkSphereSource::New();
  sphere->SetRadius(1);

  //  glyphs->SetSource(cube->GetOutput());
  glyphs->SetSource(sphere->GetOutput());

  glyphs->SetInput(targetLandmarks);
  glyphs->ScalingOff();

  vtkPolyData   *output2 = vtkPolyData::New();
  output2 = glyphs->GetOutput();

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInput(output2);


  vtkPolyDataWriter *writer2 = vtkPolyDataWriter::New();
  writer2->SetInput(cleaner->GetOutput());
  writer2->SetFileName("usingGlyphs.vtk");
  writer2->Write();



}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
