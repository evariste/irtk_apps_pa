#if (defined HAS_VTK)

#include <irtkImage.h>
//#include <nr.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

char *input_name = NULL, *output_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << "Usage: polydataclampscalars [in] [out] [lo] [hi] <options>" << endl;
  cerr << "" << endl;
  cerr << "Options: " << endl;
  cerr << "-q    : Quiet output, only print min and max values." << endl;
  cerr << "-name : Name of scalars to use." << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  double val;
  int i, count;
  int noOfPoints;
  bool quiet = false;

  // Percentiles by default
  double lo, hi;
  // values
  double minVal, maxVal;

  if (argc < 5){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;
  lo = atof(argv[1]);
  argc--;
  argv++;
  hi = atof(argv[1]);
  argc--;
  argv++;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)) {
      argc--;
      argv++;
      scalar_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData *polys = vtkPolyData::New();

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  polys = surface_reader->GetOutput();

  noOfPoints= polys->GetNumberOfPoints();

//  vtkFloatArray *scalars = vtkFloatArray::New();
//  scalars = (vtkFloatArray*) polys->GetPointData()->GetScalars();
  vtkFloatArray *scalars;

  if (scalar_name != NULL){
    scalars = (vtkFloatArray*) polys->GetPointData()->GetArray(scalar_name);
    if (scalars == NULL){
      cerr << "Cannot retrieve scalars : " << scalar_name << endl;
      exit(1);
    }
  } else {
    scalars = (vtkFloatArray*) polys->GetPointData()->GetScalars();
  }

  cerr << "Using scalars :  " << scalars->GetName() << endl;

  gsl_vector *data = gsl_vector_alloc(noOfPoints);

  count = 0;

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    gsl_vector_set(data, i, val);
  }

  gsl_sort_vector(data);

  i = (int) round( (double) lo * (noOfPoints - 1) / 100.0);
  minVal = gsl_vector_get(data, i);
  i = (int) round( (double) hi * (noOfPoints - 1) / 100.0);
  maxVal = gsl_vector_get(data, i);;

  if (quiet){
    cout << minVal << "," << maxVal << endl;
  } else {
    cout << "No of points " << noOfPoints << endl;
    cout << "Clamping to min and max : " << minVal << ", " << maxVal << endl;
  }

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    if (val < minVal)
      val = minVal;
    if (val > maxVal)
      val = maxVal;
    scalars->SetTuple1(i, val);
  }

  scalars->Modified();

  polys->GetPointData()->AddArray(scalars);
  polys->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(polys);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  return 0;
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

