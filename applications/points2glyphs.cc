#ifdef HAS_VTK

///////////////////////////////////////////////////////////////

#include <vtkFloatArray.h>
#include <vtkCubeSource.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkGlyph3D.h>
#include <irtkImage.h>

#include <vtkMath.h>


char *input_points_name = NULL;
char *output_name = NULL;


void usage(char* name){


  cout << "" << endl;
  cout << "" << endl;
   cout << name <<  " [inputPoints.vtk] [outputPolys.vtk] <-r> <-scalar>" << endl;
  cout << "" << endl;
  cout << "Place a glyph centered at each of the input points. By default the glyph is a sphere." << endl;
  cout << "default radius is 1, can assign a scalar value if wanted." << endl;


  cout << " " << endl;
  cout << " Options:" << endl;

  cout << " -r [value]       : Use the given value as the radius." << endl;
  cout << " -auto_radius     : Automatically estimate radius using minimum pairwise distance." << endl;
  cout << " -cube            : Use a cube as the glyph instead of a (default) sphere." << endl;
  cout << " " << endl;
  cout << "" << endl;

  exit(1);
}


int main(int argc, char **argv){

  bool ok;
  double p[3], q[3];
  int i, j, n;
  double d, minD = MAXFLOAT;


  if (argc < 3)
    usage(argv[0]);

  input_points_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  bool useCube = false;
  bool autoRadius = false;
  double r = 1.0;


  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-r") == 0)){
      argc--;
      argv++;
      r = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-auto_radius") == 0)){
      argc--;
      argv++;
      autoRadius = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cube") == 0)){
      argc--;
      argv++;
      useCube = true;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }


  cout << "Reading target landmarks ... "; cout.flush();
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_points_name);
  reader->Update();


  vtkPolyData* inputPoints = vtkPolyData::New();
  inputPoints = reader->GetOutput();
  cout << "done" << endl;




  n = inputPoints->GetNumberOfPoints();


  vtkFloatArray *ptIDs = vtkFloatArray::New();
  ptIDs->SetNumberOfComponents(1);
  ptIDs->SetNumberOfTuples(n);
  ptIDs->SetName("ptID");

  for (i = 0; i < n; i++){
    ptIDs->SetTuple1(i, i+1);
  }

  inputPoints->GetPointData()->AddArray(ptIDs);
  inputPoints->Update();

  for (i = 0; i < n-1; i++){
    inputPoints->GetPoint(i, p);
    for (j = i+1; j < n; j++){
      inputPoints->GetPoint(j, q);
      d = sqrt( vtkMath::Distance2BetweenPoints(p, q) );

      if (d < minD)
        minD = d;
    }
  }


  cout << "Minimum distance between points: " << minD << endl;


  vtkGlyph3D *glyphs = vtkGlyph3D::New();


  if (autoRadius){
    r = 0.1 * minD;
  }

  cout << "Setting radius to " << r << endl;

  vtkCubeSource *cube = vtkCubeSource::New();
  cube->SetBounds(-r, r, -r, r, -r, r);

  vtkSphereSource *sphere = vtkSphereSource::New();
  sphere->SetRadius(r);

  if (useCube){
    cout << "Using cubic glyph" << endl;
    glyphs->SetSource(cube->GetOutput());
  } else {
    cout << "Using spherical glyph" << endl;
    glyphs->SetSource(sphere->GetOutput());
  }

  glyphs->SetInput(inputPoints);
  glyphs->SetColorModeToColorByScalar();
  glyphs->SetScaleFactor(1.0);

  vtkPolyData   *output = vtkPolyData::New();
  output = glyphs->GetOutput();

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInput(output);
  cleaner->Update();

  vtkPolyData *output2 = vtkPolyData::New();
  output2 = cleaner->GetOutput();
  output2->GetPointData()->SetActiveScalars("ptID");
  output2->Update();


  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(output2);
  writer->SetFileName(output_name);
  writer->Write();



}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
