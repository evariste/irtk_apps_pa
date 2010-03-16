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
  cerr << " Usage: polydatasmooth [input] [output] [iterations] [relaxationFactor]  <options>" << endl;
  cerr << "" << endl;
  cerr << " Area weighted Laplacian relaxation of a surface." << endl;
  cerr << "" << endl;
  cerr << " Iteratively move mesh nodes towards the centroid of the adjacent nodes." << endl;
  cerr << " The relaxation factor (RF) determines how far each node moves towards the " << endl;
  cerr << " local centroid (0<= RF <= 1).  The new position of a node is a weighted " << endl;
  cerr << " combination of its previous position and the neighbours' centroid.  " << endl;
  cerr << " RF = 1 moves a node all the way to the centroid while RF = 0 keeps a " << endl;
  cerr << " node at its previous position." << endl;
  cerr << "" << endl;
  cerr << "Options: " << endl;
  cerr << "" << endl;
  cerr << "-track : " << endl;
  cerr << "         Track the signed distances traversed by points (out is positive)." << endl;
  cerr << "         Store the results in scalar array \"smoothingDists\"" << endl;
  cerr << "-threshold [value] " << endl;
  cerr << "         A value indicating the smoothness as the L_2 norm of H^2 over the surface." << endl;
  cerr << "         See Tosun, MedIA, 2004.  Iterations stop if the norm for the surface drops" << endl;
  cerr << "         drops below the given value." << endl;
  cerr << "" << endl;
  cerr << "" << endl;

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
  
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(smoother->GetOutput());
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

