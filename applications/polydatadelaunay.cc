#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkDelaunay3D.h>
#include <vtkCleanPolyData.h>


char *in_name  = NULL, *out_name = NULL;

void usage(char* name)
{
  cerr << "Usage: " << name << "  [input points] [output polys]\n" << endl;
  cerr << "Uses delaunay3D, specifically aimed for landmarks around a cortex .." << endl;
  cerr << "i.e. the points are arranged on the surface of a roughly spherical shape." << endl;
  cerr << "Points internal will not end up connected with the surface." << endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if (argc < 3){
    usage(argv[0]);
  }

  // Parse source and target point lists
  in_name  = argv[1];
  argc--;
  argv++;
  out_name = argv[1];
  argc--;
  argv++;

  // Read the input points.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();

  // vtkDelaunay3D returns an unstructured grid
  // (vtkUnstructuredGridAlgorithm) with tetrahedral cells.
  vtkDelaunay3D *del3d = vtkDelaunay3D::New();
  del3d->SetInput(input);

  // geometry filter is a vtkPolyDataAlgorithm which has polydata output,
  // use it to get the boundary polygons of the result of the delaunay 3D
  // algorithm.
  vtkGeometryFilter *gem = vtkGeometryFilter::New();
  gem->SetInput(  del3d->GetOutput() );
  gem->Modified();
  gem->Update();

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInput(gem->GetOutput() );
  cleaner->Update();

  // Write the output of the geometry filter.
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(out_name);
  writer->SetInput( cleaner->GetOutput() );
  writer->Update();

  // Check that no points were lost due to being internal.
  cout << "Number of input points  : " << input->GetNumberOfPoints() << endl;
  cout << "Number of output points : " << cleaner->GetOutput()->GetNumberOfPoints() << endl;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif


