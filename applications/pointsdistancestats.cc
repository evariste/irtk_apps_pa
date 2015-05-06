/*
 * pointsdistancestats.cc
 *
 *  Created on: Feb 7, 2015
 *      Author: paulaljabar
 */




#ifdef HAS_VTK

///////////////////////////////////////////////////////////////

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkMath.h>

char *input_points_name = NULL;

void usage(char* name){


  cout << "" << endl;
  cout << "" << endl;
  cout << name <<  " [inputPoints.vtk]" << endl;
  cout << "" << endl;
  cout << " Print some statistics on the distances between pairs of points." << endl;
  cout << " Minimum, mean, maximum." << endl;
  cout << "" << endl;
  cout << "Options:" << endl;
  cout << " -verbose   : Print more than just the numbers separated by commas." << endl;

  exit(1);
}


int main(int argc, char **argv){

  bool ok, verbose;
  double p[3], q[3];
  int i, j, n, count;
  double d, sumD;
  double minD = MAXFLOAT;
  double maxD = -1.0 * MAXFLOAT;

  if (argc < 2)
    usage(argv[0]);

  input_points_name = argv[1];
  argv++;
  argc--;

  verbose = false;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-verbose") == 0)){
      argc--;
      argv++;
      verbose = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }


  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_points_name);
  reader->Update();

  vtkPolyData* inputPoints = vtkPolyData::New();
  inputPoints = reader->GetOutput();

  n = inputPoints->GetNumberOfPoints();

  sumD = 0.0;
  count = 0;
  for (i = 0; i < n-1; i++){
    inputPoints->GetPoint(i, p);

    for (j = i+1; j < n; j++){
      inputPoints->GetPoint(j, q);
      d = sqrt( vtkMath::Distance2BetweenPoints(p, q) );
      sumD += d;

      count++;

      if (d < minD)
        minD = d;
      if (d > maxD)
        maxD = d;
    }

  }

  if(verbose){
    cout << "Minimum distance between points: " << minD << endl;
    cout << "Mean    distance between points: " << sumD / (float (count)) << endl;
    cout << "Maximum distance between points: " << maxD << endl;
  } else {
    cout << minD << "," << sumD / (float (count)) << "," << maxD << endl;
  }


}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
