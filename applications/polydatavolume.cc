/*
 * polydatavolume.cc
 *
 *  Created on: Feb 6, 2015
 *      Author: paulaljabar
 */


#ifdef HAS_VTK

///////////////////////////////////////////////////////////////

#include <vtkPolyDataReader.h>
#include <vtkMassProperties.h>
#include <vtkMath.h>
#include <vtkTriangleFilter.h>

char *input_poly = NULL;


void usage(char* name){


  cout << "" << endl;
  cout << "" << endl;
   cout << name <<  " [inputPolys.vtk] " << endl;
  cout << "" << endl;
  cout << " Print the volume of the input polydata surface." << endl;
  cout << " This will fail if the input does not represent a closed surface." << endl;
  cout << "" << endl;

  exit(1);
}


int main(int argc, char **argv){

  bool ok;

  if (argc < 2)
    usage(argv[0]);

  input_poly = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-option") == 0)){
      argc--;
      argv++;
      // do stuff and perhaps
//      argc--;
//      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }


  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_poly);
  reader->Update();

  vtkTriangleFilter *triFilt = vtkTriangleFilter::New();
  triFilt->SetInputData(reader->GetOutput());
  triFilt->Update();

  vtkMassProperties *mProps = vtkMassProperties::New();
  mProps->SetInputData(triFilt->GetOutput());

  cout << mProps->GetVolume() << endl;


}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif



