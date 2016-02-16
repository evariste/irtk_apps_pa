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
  cerr << "" << endl;
  cerr << " Usage: pca_visualise_mode [meanSurface.vtk] [eigenvectors.mat] [outputSurface.vtk] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Eigenvectors should be a file describing the modes of variation around" << endl;
  cerr << " the mean surface. See pcanalysis executable." << endl;
  cerr << " The output surface has the same points as the mean surface and also" << endl;
  cerr << " has scalars to indicate the extent to which the eigenmodes deviate from " << endl;
  cerr << " the mean.  Postive values are outward and negative values are inward." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -mode n    : Which mode to use, default = 0." << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int noOfPoints;
  int mode = 0;
  int i;
  double n[3], v[3];
  bool lateralComps = false;
  
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
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-lateral") == 0)){
      argc--;
      argv++;
      lateralComps = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mode") == 0)){
      argc--;
      argv++;
      mode = atoi(argv[1]);
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

  vtkPolyData *inputRead = vtkPolyData::New();
  inputRead = reader->GetOutput();



  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SetInputData(inputRead);
  normalsFilter->AutoOrientNormalsOn();
  normalsFilter->Modified();
  normalsFilter->Update();

  vtkPolyData *input = normalsFilter->GetOutput();

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
  
  double nPart[3], lPart[3];
  double compN, compL;
  
  for (i = 0; i < noOfPoints; ++i){
  	normals->GetTuple(i, n);

  	v[0] = Ev(3*i  ,mode);
  	v[1] = Ev(3*i+1,mode);
  	v[2] = Ev(3*i+2,mode);

  	(void) vtkMath::Normalize(n);

  	if (lateralComps == true){
  	  compN = vtkMath::Dot(n, v);
      nPart[0] = compN * n[0];
      nPart[1] = compN * n[1];
      nPart[2] = compN * n[2];
      
      lPart[0] = v[0] - nPart[0];
      lPart[1] = v[1] - nPart[1];
      lPart[2] = v[2] - nPart[2];

      comps->SetTuple1(i, vtkMath::Norm(lPart));
  	  
  	} else {
      comps->SetTuple1(i, vtkMath::Dot(n, v));
  	}

  }

  char buf[255];
  sprintf(buf, "mode_%02d", mode);

  comps->SetName(buf);
  input->GetPointData()->AddArray(comps);



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

