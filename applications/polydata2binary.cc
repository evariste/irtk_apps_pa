#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydata2binary [polydata]" << endl;
  exit(1);
}

char *input_name = NULL;

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  input_name  = argv[1];
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
  writer->SetFileName(input_name);
  writer->SetFileTypeToBinary();
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
