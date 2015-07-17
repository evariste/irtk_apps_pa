/*=========================================================================

  Run Philipp Batchelor's VTK filter for getting curvature data from a
  surface. The curvature values are stored as scalars on the input polydata
	set.

	Paul Aljabar, July 2010.

=========================================================================*/


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

typedef enum { K, H, kmin, kmax} curvatureType;

void usage()
{
  cerr << "Usage: polydatacurvature [in] [out] <type>" << endl;
  cerr << "" << endl;
  cerr << "Calculate curvatures for a surface. Will add a scalar float array for all" << endl;
  cerr << "vertices containing the values. The name of the scalar array will be         " << endl;
  cerr << "one of K, H, kmin, kmax (see below).        " << endl;
  cerr << "The output polydata will have copies of any scalar arrays that pre-exist        " << endl;
  cerr << "in the input polydata.        " << endl;
  cerr << "        " << endl;
  cerr << "Options:" << endl;
  cerr << "<type>  Can be one of:" << endl;
  cerr << "        " << endl;
  cerr << "        -K      Gaussian curvature" << endl;
  cerr << "        -H      Mean curvature" << endl;
  cerr << "        -kmin   Minimum principal curvature" << endl;
  cerr << "        -kmax   Maximum principal curvature" << endl;
  cerr << "        " << endl;
  cerr << "-invert : Invert the sign of mean curvature" << endl;
  cerr << "        " << endl;
  cerr << "        The default option is Gaussian." << endl;
  cerr << "        " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int i;
  bool smoothOn = true;
  bool invert = false;
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
    ok = false;
    if ((!ok) && (strcmp(argv[1], "-K") == 0)) {
      argc--;
      argv++;
      cType = K;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-H") == 0)) {
      argc--;
      argv++;
      cType = H;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-kmin") == 0)) {
      argc--;
      argv++;
      cType = kmin;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-kmax") == 0)) {
      argc--;
      argv++;
      cType = kmax;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-smoothOff") == 0)) {
      argc--;
      argv++;
      smoothOn = false;
      ok = true;
    }
    if ((!ok) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = true;
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData *input;

  // Read surface
  vtkPolyDataReader *reader1 = vtkPolyDataReader::New();
  reader1->SetFileName(input_name);
  reader1->Modified();
  reader1->Update();
  input = reader1->GetOutput();

  vtkCurvatures *curve = vtkCurvatures::New();


  if (smoothOn){
    cout << "Smoothing ..." << endl;
    vtkSmoothPolyDataFilter *smooth = vtkSmoothPolyDataFilter::New();
    smooth->SetNumberOfIterations(smoothIterations);
    smooth->SetConvergence(smoothConvergence);
    smooth->SetInput(input);
  	curve->SetInput(smooth->GetOutput());
  } else {
  	curve->SetInput(input);
  }

  switch (cType){
    case K:
      cout << "Setting curvature to Gaussian" << endl;
      curve->SetCurvatureTypeToGaussian();
      break;
    case H:
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


  vtkPolyData *output;
  vtkPolyDataReader *reader2 = vtkPolyDataReader::New();
  reader2->SetFileName(input_name);
  reader2->Modified();
  reader2->Update();
  output = reader2->GetOutput();

  vtkFloatArray *scalars = NULL;

  // Type requested is one of K, H, kmin or kmax.

  if (cType == kmin){
    scalars = static_cast<vtkFloatArray*>(curve->GetOutput()->GetPointData()->GetScalars("Minimum_Curvature"));
    scalars->SetName("kmin");
  } else if (cType == kmax) {
    scalars = static_cast<vtkFloatArray*>(curve->GetOutput()->GetPointData()->GetScalars("Maximum_Curvature"));
    scalars->SetName("kmax");
  } else if (cType == H) {
    scalars = static_cast<vtkFloatArray*>(curve->GetOutput()->GetPointData()->GetScalars("Mean_Curvature"));
    scalars->SetName("H");
  } else if (cType == K) {
    scalars = static_cast<vtkFloatArray*>(curve->GetOutput()->GetPointData()->GetScalars("Gauss_Curvature"));
    scalars->SetName("K");
  }

  output->GetPointData()->AddArray(scalars);



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
