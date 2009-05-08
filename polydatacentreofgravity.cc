#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
// Find the average distance of a set of points from a surface 
// 

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <irtkImage.h>

char *surface_name = NULL;

int main(int argc, char **argv){

  int ok;

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
      ok = True;
//     }

    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }


  cout << "Reading data ... "; cout.flush();
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(surface_name);
  reader->Update();

  vtkPolyData* surface = vtkPolyData::New();
  surface = reader->GetOutput();
  cout << "done" << endl;

  int n = surface->GetNumberOfPoints();
  double cofgx, cofgy, cofgz;
  double point[3];

  cofgx = cofgy = cofgz = 0.0;


  for (int i = 0; i < n; i++){
    surface->GetPoint (i, point);

    cofgx += point[0];
    cofgy += point[1];
    cofgz += point[2];

  }

  if (n > 0){
    cofgx /= n;
    cofgy /= n;
    cofgz /= n;
  }

  cout << "C of G : (" << cofgx << ", " << cofgy << ", " << cofgz << ")" <<  endl;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
