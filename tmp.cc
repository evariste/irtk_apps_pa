////////////////////////////////////////////////////////////////////////
// polydatasmooth with an option to stop when the L2 norm of H (extrinsic
// curvature reaches a certain level indicating a threshold on smoothness..


#if (defined HAS_VTK)

#include <vtkPolyDataSmoothCustom.h>
#include <irtkImage.h>




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
  int ok;
  int noOfIterations;
  double relaxationFactor;
  double smoothnessThreshold = -1.0f;
  int trackingOn = False;

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
  
  vtkCurvatures *curve = vtkCurvatures::New();
  curve->SetInput(smoothedInput);
  curve->SetCurvatureTypeToGaussian();
  curve->Update();
  
  
  vtkPolyData *output = vtkPolyData::New();
  output = curve->GetOutput();
  
  
  
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

