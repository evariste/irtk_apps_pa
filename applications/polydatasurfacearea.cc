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
char *output_name = NULL;

void usage()
{
  cerr << " polydatasurfacearea [input] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Print the area of the polydata file." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " " << endl;
  cerr << " -output [file] : Save a VTK polydata file with a scalar array on the vertices" << endl;
  cerr << "                  showing the area. The output is a triangulated version" << endl;
  cerr << "                  of the input. Each vertex receives the sum of 1/3 of the" << endl;
  cerr << "                  areas of triangles it is a part of." << endl;
  cerr << "                  " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  double val;
  bool ok;

  input_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
      argc--;
      argv++;
      output_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
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

  vtkPolyData *pd = vtkPolyData::New();
  pd = triFilter->GetOutput();
  pd->Update();

  vtkCellArray* facets = pd->GetPolys();
  vtkTriangle* facet = vtkTriangle::New();

  double A = 0.0;
  double v0[3], v1[3], v2[3];

  vtkFloatArray *scalarsOut = vtkFloatArray::New();

  scalarsOut->SetNumberOfComponents(1);
  scalarsOut->SetNumberOfTuples(pd->GetNumberOfPoints());
  scalarsOut->SetName("area");

  vtkIdType f, *vert=0;
  facets->InitTraversal();
  while (facets->GetNextCell(f,vert)){

    pd->GetPoint(vert[0],v0);
    pd->GetPoint(vert[1],v1);
    pd->GetPoint(vert[2],v2);

    val = double(facet->TriangleArea(v0,v1,v2));
    A += val;

    for (int i = 0; i < 3; i++)
      scalarsOut->SetTuple1(vert[i], val / 3.0);
    }

  cout << A << endl;

  if (output_name != NULL){
    pd->GetPointData()->AddArray(scalarsOut);
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetFileName(output_name);
    writer->SetInput(pd);
    writer->Write();
  }

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
