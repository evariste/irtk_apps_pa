#if (defined HAS_VTK)

#include <irtkImage.h>

#include <nr.h>
#include <time.h>

#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

char *input_name  = NULL;
char *output_name = NULL;

void usage()
{
  cerr << " " << endl;
  cerr << " Usage: polydatarandomscalars [input] [output] " << endl; // <options>" << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  long i;

  int noOfPoints;
  time_t seconds;
  long ran2Seed;
  long ran2initialSeed;

  if (argc < 3){
    usage();
  }

  // Prepare for random stuff.
  seconds = time(NULL);
  ran2Seed = seconds;
  ran2initialSeed = -1 * ran2Seed;
  (void) ran2(&ran2initialSeed);

  // Parse arguments.
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  cerr << "Input      : " << input_name << endl;
  cerr << "Output     : " << output_name << endl;

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData* input = reader->GetOutput();
  input->Update();

	noOfPoints = input->GetNumberOfPoints();

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(noOfPoints);

  for (i = 0; i < noOfPoints; ++i){
    scalars->SetTuple1(i, ran2(&ran2Seed));
  }

  scalars->SetName("Random");
  input->GetPointData()->AddArray(scalars);

  // Write the result.
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(input);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Write();

}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
