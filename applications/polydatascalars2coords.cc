#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>

char *input_name = NULL, *output_name = NULL;
char *scalar_nameX = NULL;
char *scalar_nameY = NULL;
char *scalar_nameZ = NULL;

void usage()
{
  cerr << "Usage: polydatascalars2coords [in] [out] [scalarX] [scalarY] [scalarZ] <options>" << endl;
  cerr << "" << endl;
  cerr << "Use the scalars of surface [in] as coordinates for the locations of its points." << endl;
  cerr << "After applying the scalars, write the new surface to [out]." << endl;
  cerr << "Surface in must have single component scalars with names [scalarX] [scalarY] [scalarZ]." << endl;
  cerr << "" << endl;
  cerr << "Options: " << endl;
  cerr << "-norm radius   Normalise all coordinates to radius 1." << endl;
  cerr << "-radius r      Normalise all coordinates to radius r." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int i, j;
  int noOfPoints;
  int ind;
  double p[3], len;
  int normalise = false;
  double radius = 1.0;
  
  if (argc < 6){
    usage();
  }

  // Parse image
  input_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  scalar_nameX = argv[1];
  argc--;
  argv++;
  scalar_nameY = argv[1];
  argc--;
  argv++;
  scalar_nameZ = argv[1];
  argc--;
  argv++;

  
  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-norm") == 0)){
      argc--;
      argv++;
      normalise = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-radius") == 0)){
      argc--;
      argv++;
      normalise = true;
      radius = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
     if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }


  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();



  noOfPoints= input->GetNumberOfPoints();

  vtkFloatArray *scalarsX, *scalarsY, *scalarsZ;

  scalarsX = (vtkFloatArray*) input->GetPointData()->GetArray(scalar_nameX, ind);
  if (ind == -1 || scalarsX == NULL){
    cerr << "Scalars unavailable with name " << scalar_nameX << endl;
    exit(0);
  }
  scalarsY = (vtkFloatArray*) input->GetPointData()->GetArray(scalar_nameY, ind);
  if (ind == -1 || scalarsY == NULL){
    cerr << "Scalars unavailable with name " << scalar_nameY << endl;
    exit(0);
  }
  scalarsZ = (vtkFloatArray*) input->GetPointData()->GetArray(scalar_nameZ, ind);
  if (ind == -1 || scalarsZ == NULL){
    cerr << "Scalars unavailable with name " << scalar_nameZ << endl;
    exit(0);
  }

  if (scalarsX->GetNumberOfComponents() > 1){
    cerr << "Scalars " << scalarsX->GetName() << " has more than one component." << endl;
    exit(1);
  }
  if (scalarsY->GetNumberOfComponents() > 1){
    cerr << "Scalars " << scalarsY->GetName() << " has more than one component." << endl;
    exit(1);
  }
  if (scalarsZ->GetNumberOfComponents() > 1){
    cerr << "Scalars " << scalarsZ->GetName() << " has more than one component." << endl;
    exit(1);
  }

  for (i = 0; i < noOfPoints; ++i){
    p[0] = scalarsX->GetTuple1(i);
    p[1] = scalarsY->GetTuple1(i);
    p[2] = scalarsZ->GetTuple1(i);
    
    if (normalise == true){
      len = vtkMath::Norm(p);

      for (j = 0; j < 3; ++j){
        p[j] = p[j] * radius / len;
      }
    }

    input->GetPoints()->SetPoint(i,p);
  }
  


  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(input);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  return 0;
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

