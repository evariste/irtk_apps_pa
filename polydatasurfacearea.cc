#ifdef HAS_VTK

#include <irtkImage.h>

#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkTriangleFilter.h>

#define MAXVALS 100

char *input_name = NULL;

void usage()
{
  cerr << " polydatasurfacearea [input] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  int i;
  int ok;

  input_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-XXX") == 0)){
      argc--;
      argv++;
      // Do stuff and maybe argv++ etc.
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData *inputRead = vtkPolyData::New();
  inputRead = reader->GetOutput();
  inputRead->Update();


  vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
  triFilter->SetInput(inputRead);
  triFilter->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = triFilter->GetOutput();
  input->Update();

  vtkCellArray* facets = input->GetPolys();
  vtkTriangle* facet = vtkTriangle::New();

  double A = 0.0;
  double v0[3], v1[3], v2[3];

  vtkIdType f, *vert=0;
  facets->InitTraversal();
  while (facets->GetNextCell(f,vert)){

    input->GetPoint(vert[0],v0);
    input->GetPoint(vert[1],v1);
    input->GetPoint(vert[2],v2);

    A += double(facet->TriangleArea(v0,v1,v2));

  }

  cout << A << endl;

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
