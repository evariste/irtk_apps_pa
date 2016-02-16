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
char *source_name = NULL, *lookup_name = NULL;

void usage()
{
  cerr << "Usage: polydataremesh [target] [source] [output] <options>" << endl;
  cerr << "" << endl;
  cerr << "Apply the mesh topology of [target] to the geometry of the surface represented" << endl;
  cerr << "by [source].  I.e. For each point in [target] mesh, look up the nearest location" << endl;
  cerr << "in surface [source] and apply its coordinates to the target point." << endl;
  cerr << "" << endl;
  cerr << "Options are one or more of the following:" << endl;
  cerr << "<-dofin File>      Transformation mapping locations in [target] to locations in [source]." << endl;
  cerr << "<-lookup Surface>  Surface to lookup closest points instead of [source]." << endl;
  cerr << "                   Must have same number of points as [source]." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i;
  bool ok;
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
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-lookup") == 0)) {
      argc--;
      argv++;
      lookup_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (lookup_name == NULL){
    lookup_name = source_name;
  }
  
  cout << "Using " << lookup_name << " as a replacement file." << endl;
  
  vtkPolyDataReader *lookup_reader = vtkPolyDataReader::New();
  lookup_reader->SetFileName(lookup_name);
  lookup_reader->Modified();
  lookup_reader->Update();
  vtkPolyData *lookupSurf = lookup_reader->GetOutput();

  if (lookupSurf->GetNumberOfPoints() != sourceSurf->GetNumberOfPoints()){
    cerr << "Error : lookup surface must have same number of points as source surface. " << endl;
    exit(1);
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

  int closestID;
  
  for (i = 0; i < noOfPoints; i++) {
    tgtPoints->GetPoint(i, coord);
    transformation->Transform(coord[0], coord[1], coord[2]);
    closestID = source_locator->FindClosestPoint(coord);
    lookupSurf->GetPoint(closestID, coord);
    tgtPoints->SetPoint(i, coord);
  }

  targetSurf->SetPoints(tgtPoints);
  targetSurf->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(targetSurf);
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
