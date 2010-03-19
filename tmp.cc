////////////////////////////////////////////////////////////////////////
// polydatasmooth with an option to stop when the L2 norm of H (extrinsic
// curvature reaches a certain level indicating a threshold on smoothness..


#if (defined HAS_VTK)

#include <vtkFillHolesFilter.h>
#include <vtkPolyDataSmoothCustom.h>
#include <irtkImage.h>
#include <nr.h>




char *input_name = NULL;
char *output_name = NULL;

void usage()
{
  cerr << "" << endl;
  cerr << " Usage:" << endl;
  exit(1);
}


int main(int argc, char **argv)
{
  int ok, i, j, ind, noOfPoints;
  int noOfIterations;
  double relaxationFactor;
  double smoothnessThreshold = -1.0f;
  int trackingOn = False;
  int lo, hi;
  double val, minVal, maxVal;
  double holeSize = 50;

  lo = 2;
  hi = 98;

  if (argc < 4){
    usage();
  }

  // Parse arguments.
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;
  noOfIterations = atoi(argv[1]);
  argv++;
  argc--;
  relaxationFactor = atof(argv[1]);
  argv++;
  argc--;

  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-track") == 0)){
      argc--;
      argv++;
      trackingOn = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;
      smoothnessThreshold = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-holesize") == 0)){
      argc--;
      argv++;
      holeSize = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }    if ((ok == False) && (strcmp(argv[1], "-lohi") == 0)){
      argc--;
      argv++;
      lo = atoi(argv[1]);
      argc--;
      argv++;
      hi = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
     if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Input        : " << input_name << endl;
  cerr << "Output       : " << output_name << endl;
  cerr << "Iterations   : " << noOfIterations << endl;
  cerr << "Relax factor : " << relaxationFactor << endl;

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);

  vtkPolyData* input = reader->GetOutput();
  input->Update();

  vtkPolyDataSmoothCustom *smoother = new vtkPolyDataSmoothCustom;
  smoother->SetInput(input);
  smoother->SetNoOfIterations(noOfIterations);
  smoother->SetRelaxationFactor(relaxationFactor);
  smoother->SetSmoothnessThreshold(smoothnessThreshold);
  smoother->SetTrackingOn(trackingOn);
  smoother->Run();

  vtkPolyData *smoothedInput = smoother->GetOutput();
  smoothedInput->Update();

  vtkCurvatures *curve = vtkCurvatures::New();
  curve->SetInput(smoothedInput);
  curve->SetCurvatureTypeToGaussian();
  curve->Update();


  vtkPolyData *curveOutput = vtkPolyData::New();
  curveOutput = curve->GetOutput();
  curveOutput->Update();


  vtkFloatArray *scalars;
  scalars = (vtkFloatArray*) curveOutput->GetPointData()->GetArray("Gauss_Curvature", ind);

  if (ind == -1 || scalars == NULL){
    cerr << "Curvatures unavailable" << endl;
    exit(0);
  }

  noOfPoints = curveOutput->GetNumberOfPoints();

  float *data = new float[1 + noOfPoints];

  for (i = 0; i < noOfPoints; ++i){
    data[i+1] = scalars->GetTuple1(i);
  }

  sort(noOfPoints, data);


  i = 1 + (int) round( (double) lo * (noOfPoints - 1) / 100.0);
  minVal = data[i];
  i = 1 + (int) round( (double) hi * (noOfPoints - 1) / 100.0);
  maxVal = data[i];

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    if (val < minVal || val > maxVal){
      scalars->SetTuple1(i, 1);
    } else {
      scalars->SetTuple1(i, 0);
    }
  }

  //////////////////////////////////////
  // Delete points:

  vtkPoints *oldPoints = vtkPoints::New();
  vtkPoints *newPoints = vtkPoints::New();

  oldPoints = input->GetPoints();

  int *newPtIds, *deletedPt;
  int newPtCount;
  double p[3];

  newPtIds = new int[noOfPoints];
  deletedPt = new int[noOfPoints];

  for (i = 0; i < noOfPoints; ++i){
    deletedPt[i] = 0;
  }

  newPtCount = 0;
  for (i = 0; i < noOfPoints; ++i){

    if (scalars->GetTuple1(i) > 0){
      // Outside curvature limits. Ignore this point.
      continue;
    }
    ++newPtCount;
  }

  newPoints->SetNumberOfPoints(newPtCount);

  newPtCount = 0;
  for (i = 0; i < noOfPoints; ++i){
    if (scalars->GetTuple1(i) > 0){
      deletedPt[i] = 1;
    } else {
      // Can keep the point.
      oldPoints->GetPoint(i, p);
      newPoints->InsertPoint(newPtCount, p);
      newPtIds[i] = newPtCount;
      ++newPtCount;
    }
  }

  cerr << "Points in original surface: " << noOfPoints << endl;
  cerr << "Points copied over        : " << newPtCount << endl;

  ////////////////////////////////
  // Delete faces.

  vtkCellArray* oldFaces = vtkCellArray::New();
  vtkCellArray* newFaces = vtkCellArray::New();

  oldFaces = input->GetPolys();

  int noOfFaces, newFaceCount, deleteFace, newId;
  vtkIdType npts = 0;
  vtkIdType *oldPtId;

  noOfFaces = oldFaces->GetNumberOfCells();

  newFaces->Initialize();
  newFaces->Allocate(noOfFaces, 0);
  newFaces->Reset();

  oldFaces->InitTraversal();
  newFaces->InitTraversal();
  newFaceCount = 0;

  for (i = 0; i < noOfFaces; ++i){
    oldFaces->GetNextCell(npts, oldPtId);

    deleteFace = 0;
    for (j = 0; j < npts; ++j){
      if(deletedPt[oldPtId[j]] == 1){
        deleteFace = 1;
      }
    }

    if (deleteFace == 1)
      continue;

    newFaces->InsertNextCell(npts);
    for (j = 0; j < npts; ++j){
      newId = newPtIds[oldPtId[j]];
      newFaces->InsertCellPoint(newId);
    }
    newFaces->UpdateCellCount(npts);
    ++newFaceCount;
  }
  newFaces->Squeeze();

  cerr << "Faces in original surface : " << noOfFaces << endl;
  cerr << "Faces copied over         : " << newFaceCount << endl;


  /////////////////////////////////////////////////////


  vtkPolyData* polyWithHoles = vtkPolyData::New();
  polyWithHoles->SetPoints(newPoints);
  polyWithHoles->SetPolys(newFaces);
  polyWithHoles->Update();

  vtkFillHolesFilter *holeFiller = vtkFillHolesFilter::New();
  holeFiller->SetInput(polyWithHoles);
  holeFiller->SetHoleSize(holeSize);



  vtkPolyData* output = vtkPolyData::New();
  output = holeFiller->GetOutput();
  output->Update();


  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(output);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Write();

  return 0;
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

