#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydataconvert [polydataIN] [polydataOUT] [-a | -b]" << endl;
  cerr << "Where :" << endl;
  cerr << "  -a  = convert to ascii." << endl;
  cerr << "  -b  = convert to binary." << endl;
  exit(1);
}

char *input_name = NULL;
char *output_name = NULL;
char *option_name = NULL;

int main(int argc, char **argv)
{
  if (argc < 3){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;

  output_name  = argv[1];
  argc--;
  argv++;

  option_name  = argv[1];
  argc--;
  argv++;

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(output_name);

  if (strcmp(option_name, "-a") == 0){
    writer->SetFileTypeToASCII();
  }
  if (strcmp(option_name, "-b") == 0){
    writer->SetFileTypeToBinary();
  }

  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
