#ifdef HAS_VTK

#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>

char *input_name, *output_name;

void usage()
{
  cerr << "Usage: polydatadecimate [polydataIn] [polydataOut] <-reduction val> <-smooth>" << endl;
  cerr << "-reduction val     Value between 0 and 1.  New surface should have about" << endl;
  cerr << "                   val x number of points in input surface (default = 0.5)." << endl;
  cerr << "-smooth            Smooth the surface after decimation. " << endl;
  cerr << "-ascii             Write output as ascii, default = binary." << endl;
  exit(1); 
}

int main(int argc, char **argv)
{
  int ok;
  double reduction = 0.5;
  vtkSmoothPolyDataFilter *smooth = NULL;  

  int fileType = VTK_BINARY;


  if (argc < 3){
    usage();
  }
  
  // Parse parameters
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Modified();
  reader->Update();
  vtkPolyData *surface = reader->GetOutput();

  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((!ok) && (strcmp(argv[1], "-smooth") == 0)) {
      argc--;
      argv++;
      smooth = vtkSmoothPolyDataFilter::New();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-reduction") == 0)) {
      argc--;
      argv++;
      reduction = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-ascii") == 0)) {
      argc--;
      argv++;
      fileType = VTK_ASCII;
      ok = True;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkDecimatePro *decimate = vtkDecimatePro::New();
  decimate->SetTargetReduction(reduction);
  decimate->SetInput(surface);

  if (smooth != NULL) {
    cout << "Smoothing ... \n";
    smooth->SetInput(decimate->GetOutput());
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(smooth->GetOutput());
    writer->SetFileName(output_name);
    writer->SetFileType(fileType);
    writer->Write();
  } else {
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(decimate->GetOutput());
    writer->SetFileName(output_name);
    writer->SetFileType(fileType);
    writer->Write();
  }

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}

#endif

