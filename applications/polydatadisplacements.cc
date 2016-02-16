#ifdef HAS_VTK

#include <irtkImage.h>

#include <irtkTransformation.h>

#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

// Default filenames
char *input_name = NULL, *output_name = NULL, *dof_name  = NULL;

void usage()
{
  cerr << "Usage: polydatadisplacements [surfaceIn] [surfaceOut] [dof] <options>\n" << endl;
  cerr << "  " << endl;
  cerr << "  Find the displacement magnitudes of [dof] at the points of [surfaceIn]" << endl;
  cerr << "  Apply them as a scalar 'distances' and write the result to [surfaceOut]" << endl;
  cerr << "  " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i;
  bool ok;
  int noOfPoints;

  irtkTransformation *transformation = NULL;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-option") == 0)) {
      argc--;
      argv++;
//Do stuff
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  // Read transformation
  transformation = irtkTransformation::New(dof_name);


  noOfPoints = surface->GetNumberOfPoints();

  vtkFloatArray *dists = vtkFloatArray::New();
  dists->SetNumberOfComponents(1);
  dists->SetNumberOfTuples(noOfPoints);

  double ptA[3], ptB[3];

  for (i = 0; i < noOfPoints; ++i) {
    (surface->GetPoints())->GetPoint(i, ptA);
    (surface->GetPoints())->GetPoint(i, ptB);
    transformation->Transform(ptB[0], ptB[1], ptB[2]);

    dists->SetTuple1(i, vtkMath::Distance2BetweenPoints(ptA, ptB));
  }

  dists->SetName("distance");
  surface->GetPointData()->AddArray(dists);
  surface->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(surface);
  writer->SetFileName(output_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
