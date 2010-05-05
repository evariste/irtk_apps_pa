#if (defined HAS_VTK)

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkCurvatures.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDecimatePro.h>

char *input_name = NULL, *output_name = NULL;

typedef enum { K, H, H2, kmin, kmax, Kabs, K2 } curvatureType;

void usage()
{
  cerr << "Usage: polydatacurvature [in] [out] <type>" << endl;
  cerr << "" << endl;
  cerr << "Calculate curvatures for a surface." << endl;
  cerr << "Options:" << endl;
  cerr << "<type>  Can be one of {-K -H -H2 -kmin -kmax -Kabs -K2} representing the " << endl;
  cerr << "        Gaussian curvature, mean curvature, the two principal curvatures," << endl;
  cerr << "        the absolute Gaussian and the squared Gaussian curvature." << endl;
  cerr << "        The default option is Gaussian." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i;
  int smoothOn = True;
  int decimateOn = False;
  double decimateTarget = 0;
  int cleanOn = True;
  int invert = False;
  int invertFlag;
  int smoothIterations = 200;
  double smoothConvergence = 0.001;
  int noOfPoints;
  double val;

  if (argc < 3){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  curvatureType cType = K;

  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((!ok) && (strcmp(argv[1], "-K") == 0)) {
      argc--;
      argv++;
      cType = K;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-H") == 0)) {
      argc--;
      argv++;
      cType = H;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-H2") == 0)) {
      argc--;
      argv++;
      cType = H2;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-kmin") == 0)) {
      argc--;
      argv++;
      cType = kmin;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-kmax") == 0)) {
      argc--;
      argv++;
      cType = kmax;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-Kabs") == 0)) {
      argc--;
      argv++;
      cType = Kabs;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-K2") == 0)) {
      argc--;
      argv++;
      cType = K2;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-dec") == 0)) {
      argc--;
      argv++;
      decimateOn = True;
      decimateTarget = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-smoothOff") == 0)) {
      argc--;
      argv++;
      smoothOn = False;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-cleanOff") == 0)) {
      argc--;
      argv++;
      cleanOn = False;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = True;
      ok = True;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData *surface;

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  surface = surface_reader->GetOutput();

  if (1 == 0){ // Removed, cleaning makes bad things happen to the surface //  (cleanOn){
    cout << "Cleaning ..." << endl;
    vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
    cleaner->SetInput(surface);
    cleaner->SetTolerance(0.005);
    surface = cleaner->GetOutput();
  }

  if (decimateOn){
    cout << "Decimating ..." << endl;
    vtkDecimatePro *decimate = vtkDecimatePro::New();
    decimate->SetInput(surface);
    decimate->SetTargetReduction(decimateTarget);
    surface = decimate->GetOutput();
  }

  if (smoothOn){
    cout << "Smoothing ..." << endl;
    vtkSmoothPolyDataFilter *smooth = vtkSmoothPolyDataFilter::New();
    smooth->SetNumberOfIterations(smoothIterations);
    smooth->SetConvergence(smoothConvergence);
    smooth->SetInput(surface);
    surface = smooth->GetOutput();
  }

  vtkCurvatures *curve = vtkCurvatures::New();
  curve->SetInput(surface);

  switch (cType){
    case K:
    case Kabs:
    case K2:
      cout << "Setting curvature to Gaussian" << endl;
      curve->SetCurvatureTypeToGaussian();
      break;
    case H:
    case H2:
      cout << "Setting curvature to mean" << endl;
      curve->SetCurvatureTypeToMean();
      break;
    case kmin:
      cout << "Setting curvature to kmin" << endl;
      curve->SetCurvatureTypeToMinimum();
      break;
    case kmax:
      cout << "Setting curvature to kmax" << endl;
      curve->SetCurvatureTypeToMaximum();
      break;
    default:
      cerr << "Invalid curvature type." << endl;
      exit(1);
  }

  if (invert){
    cout << "Inverting mean curvature flag." << endl;
    invertFlag = curve->GetInvertMeanCurvature();
    curve->SetInvertMeanCurvature(1 - invertFlag);
  }

  // Retrieve the computed curvature values.
  curve->Update();

//  vtkPolyData *output1 = vtkPolyData::New();
//  output1 = curve->GetOutput();
//  output1->Update();
//
//  vtkCurvatures *curve2 = vtkCurvatures::New();
//  curve2->SetInput(surface);
//  curve2->SetCurvatureTypeToMean();
//  curve2->Update();
//  vtkPolyData *output2 = vtkPolyData::New();
//  output2 = curve2->GetOutput();
//  output2->Update();
//
//  vtkFloatArray *scalars = vtkFloatArray::New();
//  output1->GetPointData()->SetActiveScalars("Gauss_Curvature");
//  scalars = (vtkFloatArray*) output1->GetPointData()->GetScalars();
//
//  cerr << "bla" << endl;
//  cerr << scalars->GetName() << endl;
//  cerr << "-------------" << endl;
//  for (i = 0; i < output1->GetPointData()->GetNumberOfArrays(); ++i){
//    cerr << output1->GetPointData()->GetArray(i)->GetName() << endl;
//  }
//  cerr << "-------------" << endl;
//
//  output2->GetPointData()->AddArray(scalars);
//  output2->Update();

  vtkPolyData *output = vtkPolyData::New();
  output = curve->GetOutput();

  vtkFloatArray *scalars = vtkFloatArray::New();

  // Any more work to be done?
  if (cType == K2 || cType == Kabs){

    scalars = (vtkFloatArray *) output->GetPointData()->GetScalars("Gauss_Curvature");
    noOfPoints = scalars->GetNumberOfTuples();

    if (cType == K2){
      scalars->SetName("K2");
      cout << "Squaring  ..." << endl;
      for (i = 0; i < noOfPoints; ++i){
        val = scalars->GetTuple1(i);
        scalars->SetTuple1(i, val * val);
      }
    } else if (cType == Kabs){
      scalars->SetName("K_abs");
      cout << "Finding absolute value ..." << endl;
      for (i = 0; i < noOfPoints; ++i){
        val = scalars->GetTuple1(i);
        scalars->SetTuple1(i, fabs(val));
      }
    }
    output->GetPointData()->AddArray(scalars);
    output->Update();
    output->GetPointData()->RemoveArray("Gauss_Curvature");
    output->Update();

    if (cType == K2){
      output->GetPointData()->SetActiveScalars("K2");
    } else {
      output->GetPointData()->SetActiveScalars("K_abs");
    }


  } else if (cType == H2){

    scalars = (vtkFloatArray *) output->GetPointData()->GetScalars();
    noOfPoints = scalars->GetNumberOfTuples();

    scalars->SetName("H2");
    cout << "Squaring  ..." << endl;

    for (i = 0; i < noOfPoints; ++i){
      val = scalars->GetTuple1(i);
      scalars->SetTuple1(i, val * val);
    }

    output->GetPointData()->AddArray(scalars);
    output->Update();
    output->GetPointData()->RemoveArray("Mean_Curvature");
    output->Update();
    output->GetPointData()->SetActiveScalars("H2");
  }

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();

  writer->SetInput(output);
//  writer->SetInput(output2);

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
