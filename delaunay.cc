#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkDelaunay3D.h>

char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: delaunay [input points] [output polys]\n" << endl;
  cerr << "Uses delaunay3D, specifically aimed for landmarks around a cortex .." << endl;
  cerr << "i.e. the points are arranged on the surface of a roughly spherical shape." << endl;
  cerr << "Points internal will not end up connected with the surface." << endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if (argc < 3){
    usage();
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

  // Write the output of the geometry filter.
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(out_name);
  writer->SetInput( gem->GetOutput() );
  writer->Update();

  // Check that no points were lost due to being internal.
  cout << "Number of input points  : " << input->GetNumberOfPoints() << endl;
  cout << "Number of output points : " << gem->GetOutput()->GetNumberOfPoints() << endl;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

//   int nPoints = input->GetNumberOfPoints();

//   cout << "nPoints " << nPoints << endl;

//   double meanZ = 0.0;
//   double xyz[3];

//   for (int i = 0; i < nPoints; i++){
//     input->GetPoint(i, xyz);
//     meanZ += xyz[2];
//   }

//   meanZ /= nPoints;
//   cout << "meanZ " << meanZ << endl;

//   // to contain float by default.
//   vtkPoints *pointsUp  = vtkPoints::New();
//   vtkPoints *pointsDown = vtkPoints::New();

//   for (int i = 0; i < nPoints; i++){
//     input->GetPoint(i, xyz);
//     if (xyz[2] < meanZ){
//       pointsDown->InsertNextPoint(xyz);
//     } else {
//       pointsUp->InsertNextPoint(xyz);
//     }
//   }

//   cout << "Points up  : " << pointsUp->GetNumberOfPoints() << endl;
//   cout << "Points down: " << pointsDown->GetNumberOfPoints() << endl;

//   vtkPolyData *polysUp   = vtkPolyData::New();
//   vtkPolyData *polysDown = vtkPolyData::New();

//   polysUp->SetPoints(pointsUp);

//   Delaunay2D , gives output of type polydata
//   vtkDelaunay2D *delaunay = vtkDelaunay2D::New();
//   delaunay->SetInput(polysUp);
//    delaunay->SetTolerance(0.01);
//    delaunay->SetAlpha(2);
//    delaunay->BoundingTriangulationOff();

//   polysUp = delaunay->GetOutput();
//   polysUp->Update();

//   MarkBoundary(polysUp);

//   int nPointsUp = polysUp->GetNumberOfPoints();

//   for (int i = 0; i < nPointsUp; i++){
//     if (*polysUp->GetPointData()->GetScalars()->GetTuple (id) == 0){
//       // A boundary point.
//       polysUp->GetPoint(i, xyz);
//       // Include with the lower point set.
//       pointsDown->InsertNextPoint(xyz);
//     }
//   }
//   polysDown->SetPoints(pointsDown);

//   delaunay->SetInput(polysDown);
//   polysDown = delaunay->GetOutput();
//   polysDown->Update();

//   vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//   writer->SetFileName(out_name);
//   writer->SetInput(polysUp);
//   writer->Update();





//   vtkUnsignedCharArray *charArray = vtkUnsignedCharArray::New();
//   charArray = ugrid->GetCellTypesArray();

//   unsigned char *temp;
//   // From vtkAbstractArray
//   int count = charArray->GetNumberOfTuples();
//   cout << "Count = " << count << endl;
//   int n;

//   for (int i = 0; i < count; i++){
//     charArray->GetTupleValue(i, temp);
//     cout << temp << endl;
//   }

//   vtkCell *cell = ugrid->GetCell(0);
//   cout << cell->GetClassName() << endl << endl;

//   cout << *charArray;

//   vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
//   writer->SetFileName(out_name);
//   writer->SetInput(ugrid);
//   writer->Update();

  //  vtkUnstructuredGridToPolyDataFilter *gem = vtkUnstructuredGridToPolyDataFilter::New();
//   vtkGeometryFilter *gem = vtkGeometryFilter::New();
//   gem->SetInput(  delaunay->GetOutput() );
//   gem->Modified();
//   gem->Update();



