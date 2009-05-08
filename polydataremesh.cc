#ifdef HAS_VTK

#include <irtkImage.h>

#include <irtkTransformation.h>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <irtkLocator.h>

// Default filenames
char *target_name = NULL, *output_name = NULL, *dof_name  = NULL;
char *source_name = NULL;

void usage()
{
  cerr << "Usage: polydataremesh [target] [source] [output] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-dofin file>      Transformation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, ok;
  irtkTransformation *transformation = NULL;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse image
  target_name  = argv[1];
  argc--;
  argv++;
  source_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read surface
  vtkPolyDataReader *target_reader = vtkPolyDataReader::New();
  target_reader->SetFileName(target_name);
  target_reader->Modified();
  target_reader->Update();
  vtkPolyData *targetSurf = target_reader->GetOutput();

  vtkPolyDataReader *source_reader = vtkPolyDataReader::New();
  source_reader->SetFileName(source_name);
  source_reader->Modified();
  source_reader->Update();
  vtkPolyData *sourceSurf = source_reader->GetOutput();

  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (dof_name != NULL) {
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }

  double coord[3];
  // vtkKDTreePointLocator
  int locatorType = 2;

  int noOfPoints;

  noOfPoints = targetSurf->GetNumberOfPoints();
  // Create locator
  irtkLocator *source_locator = new irtkLocator;
  source_locator->SelectLocatorType(locatorType);
  source_locator->SetDataSet(sourceSurf);

  vtkPoints *tgtPoints = vtkPoints::New();
  tgtPoints = targetSurf->GetPoints();

  for (i = 0; i < noOfPoints; i++) {
    tgtPoints->GetPoint(i, coord);
    transformation->Transform(coord[0], coord[1], coord[2]);
    (void) source_locator->FindClosestPoint(coord);
    tgtPoints->SetPoint(i, coord);
  }

  targetSurf->SetPoints(tgtPoints);
  targetSurf->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(targetSurf);
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
