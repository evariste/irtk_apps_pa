#ifdef HAS_VTK

#include <irtkImage.h>

#include <irtkTransformation.h>

#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
//#include <irtkLocator.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>
#include <vtkTriangleFilter.h>
#include <vtkMath.h>

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
  int i, j, ok;
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
  vtkPolyData *sourceRead = source_reader->GetOutput();
  
  
  vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
  triFilter->SetInput(sourceRead);
  triFilter->Update();
  vtkPolyData *sourceSurf = triFilter->GetOutput();
    

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

  double tgtPt[3], val;

  int noOfPoints;
  int ptID;

  noOfPoints = targetSurf->GetNumberOfPoints();
  // Create locator
  vtkCellLocator *source_locator = vtkCellLocator::New();
  source_locator->SetDataSet(sourceSurf);
  source_locator->BuildLocator();
  

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
    
    vtkGenericCell *cell = vtkGenericCell::New();
    int cellId, subId;
    double dist2;
    double closestPoint[3];
    vtkIdList *ptIDs;
    double srcPtA[3], srcPtB[3], srcPtC[3];
    double crossA[3], crossB[3], crossC[3];
    double tgt2srcPtA[3], tgt2srcPtB[3], tgt2srcPtC[3];
    int srcIdA, srcIdB, srcIdC;
    double wA, wB, wC, wSum;
    double valA, valB, valC;
    
    
    for (i = 0; i < noOfPoints; i++) {
      tgtPoints->GetPoint(i, tgtPt);
      transformation->Transform(tgtPt[0], tgtPt[1], tgtPt[2]);
      
      source_locator->FindClosestPoint(tgtPt, &closestPoint[0], 
          cell, (vtkIdType&) cellId, subId, dist2);
      
      if (cell->GetNumberOfPoints() != 3){
        cerr << "Error : filtered source should only have triangles." << endl;
        exit(1);
      }
      
      // Interpolation on a triangular face.
      
      ptIDs = cell->GetPointIds();
      
      srcIdA = ptIDs->GetId(0);
      srcIdB = ptIDs->GetId(1);
      srcIdC = ptIDs->GetId(2);
      
      sourceSurf->GetPoint(srcIdA, srcPtA);
      sourceSurf->GetPoint(srcIdB, srcPtB);
      sourceSurf->GetPoint(srcIdC, srcPtC);
      
      for (j = 0; j < 3; ++j){
        tgt2srcPtA[j] = srcPtA[j] - tgtPt[j];
        tgt2srcPtB[j] = srcPtB[j] - tgtPt[j];
        tgt2srcPtC[j] = srcPtC[j] - tgtPt[j];
      }
      
      vtkMath::Cross(tgt2srcPtA, tgt2srcPtB, crossC);
      vtkMath::Cross(tgt2srcPtB, tgt2srcPtC, crossA);
      vtkMath::Cross(tgt2srcPtC, tgt2srcPtA, crossB);
      
      wA = vtkMath::Norm(crossA);
      wB = vtkMath::Norm(crossB);
      wC = vtkMath::Norm(crossC);
      wSum = wA + wB + wC;
      
      valA = scalarsIn->GetTuple1(srcIdA);
      valB = scalarsIn->GetTuple1(srcIdB);
      valC = scalarsIn->GetTuple1(srcIdC);

      val = (wA * valA + wB * valB + wC * valC) / wSum;
      
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
