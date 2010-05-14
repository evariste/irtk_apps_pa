#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <vtkPolyDataNormals.h>

char *input_name = NULL, *output_name = NULL;
char *eigenvector_name = NULL;
char *scalar_nameX = NULL;
char *scalar_nameY = NULL;
char *scalar_nameZ = NULL;

void usage()
{
  cerr << "Usage: XXXX [meanSurf] [eigenvectors] [outsurf]" << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int ok, noOfPoints;
  int mode = 0;
  int i;
  double n[3], v[3];
  
  if (argc < 4){
    usage();
  }

  // Parse image
  input_name = argv[1];
  argc--;
  argv++;
  eigenvector_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  
  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-option") == 0)){
      argc--;
      argv++;

      // maybe argc-- etc.
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-mode") == 0)){
      argc--;
      argv++;
      mode = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
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

  vtkPolyData *inputRead = vtkPolyData::New();
  inputRead = reader->GetOutput();
  inputRead->Update();



  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SetInput(inputRead);
  normalsFilter->AutoOrientNormalsOn();
  normalsFilter->Modified();
  normalsFilter->Update();

  vtkPolyData *input = normalsFilter->GetOutput();
  input->Update();

  vtkFloatArray *normals = vtkFloatArray::New();
  normals = (vtkFloatArray*) input->GetPointData()->GetNormals();

  // Read eigenvectors
  irtkMatrix Ev;
  Ev.Read(eigenvector_name);


  noOfPoints= input->GetNumberOfPoints();

  // Components of Eigenmode along surface normal.
  vtkFloatArray *comps = vtkFloatArray::New();
  comps->SetNumberOfComponents(1);
  comps->SetNumberOfTuples(noOfPoints);

  for (i = 0; i < noOfPoints; ++i){
  	normals->GetTuple(i, n);

  	v[0] = Ev(3*i  ,mode);
  	v[1] = Ev(3*i+1,mode);
  	v[2] = Ev(3*i+2,mode);

  	(void) vtkMath::Normalize(n);
  	comps->SetTuple1(i, vtkMath::Dot(n, v));

  }

  char buf[255];
  sprintf(buf, "mode_%02d.vtk", mode);

  comps->SetName(buf);
  input->GetPointData()->AddArray(comps);

  input->Update();


  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(input);
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
