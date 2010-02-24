#ifdef HAS_VTK

#include <irtkImage.h>

#include <irtkTransformation.h>

#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <irtkLocator.h>

// Default filenames
char *target_name = NULL, *output_name = NULL, *dof_name  = NULL;
char *source_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << "Usage: polydatatransformscalars [target] [source] [output] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-dofin file>      Transformation" << endl;
  cerr << "<-scalar name>     Name of scalar array to be transformed." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, ok;
  int matching = False;

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
    if ((ok == False) && (strcmp(argv[1], "-scalar") == 0)) {
      argc--;
      argv++;
      scalar_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-matching") == 0)) {
      argc--;
      argv++;
      matching = True;
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

  vtkFloatArray *scalarsIn;

  if (scalar_name != NULL){
    scalarsIn = (vtkFloatArray*) sourceSurf->GetPointData()->GetArray(scalar_name);
    if (scalarsIn == NULL){
      cerr << "Cannot retrieve scalars : " << scalar_name << endl;
      exit(1);
    }
  } else {
    scalarsIn = (vtkFloatArray*) sourceSurf->GetPointData()->GetScalars();
    cerr << "Using scalars :  " << scalarsIn->GetName() << endl;
    scalar_name = scalarsIn->GetName();
  }

  double coord[3], val;
  // vtkKDTreePointLocator
  int locatorType = 2;

  int noOfPoints;
  int ptID;

  noOfPoints = targetSurf->GetNumberOfPoints();
  // Create locator
  irtkLocator *source_locator = new irtkLocator;
  source_locator->SelectLocatorType(locatorType);
  source_locator->SetDataSet(sourceSurf);

  vtkPoints *tgtPoints = vtkPoints::New();
  tgtPoints = targetSurf->GetPoints();

  vtkFloatArray *scalarsOut = vtkFloatArray::New();
  scalarsOut->SetNumberOfComponents(1);
  scalarsOut->SetNumberOfTuples(noOfPoints);

  if (matching == True){
    if (noOfPoints != sourceSurf->GetNumberOfPoints()){
      cerr << "Matching flag chosen but numbers of points unequal." << endl;
      exit(1);
    }
    for (i = 0; i < noOfPoints; i++) {
      val = scalarsIn->GetTuple1(i);
      scalarsOut->SetTuple1(i, val);
    }

  } else {
    for (i = 0; i < noOfPoints; i++) {
      tgtPoints->GetPoint(i, coord);
      transformation->Transform(coord[0], coord[1], coord[2]);
      ptID = source_locator->FindClosestPoint(coord);
      val = scalarsIn->GetTuple1(ptID);
      scalarsOut->SetTuple1(i, val);
    }
  }

  scalarsOut->SetName(scalar_name);


  targetSurf->GetPointData()->AddArray(scalarsOut);
  targetSurf->GetPointData()->SetActiveScalars(scalar_name);
  targetSurf->Modified();
  targetSurf->Update();

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
