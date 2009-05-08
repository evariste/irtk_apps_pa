#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>


char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << "polydatalcc [input] [output] <-options>"<<endl;
  cerr << "Keep only largest connected component of the input polydata." << endl;
  cerr << "Options:" << endl;
  exit(1);
}

int main(int argc, char **argv ){

  if (argc <3 ){
    usage();
  }

  int i, j, noOfPoints, regionCount, ok;
 
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-XX") == 0)){
      argc--;
      argv++;
      // do stuff and maybe increment argv
      ok  = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData* inputPoly;
  vtkPolyDataReader* reader;
  vtkPolyDataConnectivityFilter* connFilter;
  vtkCleanPolyData* cleaner;
  
  cout << "Reading file " << input_name << endl;

  inputPoly = vtkPolyData::New();
  reader    = vtkPolyDataReader::New();
  connFilter = vtkPolyDataConnectivityFilter::New();
  cleaner   = vtkCleanPolyData::New();
 
  connFilter->SetExtractionModeToLargestRegion();

  reader->SetFileName(input_name);
  reader->Update();

  inputPoly = reader->GetOutput();
  inputPoly->Update();

  connFilter->SetInput(inputPoly);
  connFilter->Update();


  cleaner->SetInput(connFilter->GetOutput());
  cleaner->Update();

  cout << "Output file " << output_name << endl;

  vtkPolyDataWriter *pd_writer = vtkPolyDataWriter::New();
  pd_writer->SetFileName(output_name);
  pd_writer->SetInput(cleaner->GetOutput());
  pd_writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
