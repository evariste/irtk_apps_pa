#if (defined HAS_VTK)

#include <irtkImage.h>
#include <irtkHistogram.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

char *input_name = NULL, *output_name = NULL;


void usage()
{
  cerr << "Usage: polydatascalarhist [in] <options>" << endl;
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "-bins" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  int printAll = False;
  double loHi[2];
  int loHiSet = False;
  double range;
  int nBins = 100;
  double binWidth;
  double val;
  int i;
  int noOfPoints;

  if (argc < 2){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  

  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-f") == 0)){
      argc--;
      argv++;
      output_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-all") == 0)){
      argc--;
      argv++;
      printAll = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bins") == 0)){
      argc--;
      argv++;
      nBins = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-limits") == 0)){
      argc--;
      argv++;
      loHi[0] = atof(argv[1]);
      argc--;
      argv++;
      loHi[1] = atof(argv[1]);
      argc--;
      argv++;
      loHiSet = True;
      ok = True;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkHistogram_1D hist;
  hist.Reset();

  vtkPolyData *surface = vtkPolyData::New();

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  surface = surface_reader->GetOutput();
  
  noOfPoints= surface->GetNumberOfPoints();
  cout << "No of points " << noOfPoints << endl;

  vtkFloatArray *scalars = vtkFloatArray::New();

  scalars = (vtkFloatArray*) surface->GetPointData()->GetScalars();

  if (loHiSet == False){
    surface->GetPointData()->GetScalars()->GetRange(loHi);
  }

  range = loHi[1] - loHi[0];
  binWidth = range / ((double) nBins);// - 1.0);

  hist.PutMin(loHi[0]);// - 0.5 * binWidth);
  hist.PutMax(loHi[1]);// + 0.5 * binWidth);

  cout << "lo,hi:         " << loHi[0] << "     " << loHi[1] << endl;


  hist.PutNumberOfBins(nBins);

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    hist.AddSample(val);
  }

  // Write output
  if (output_name != NULL){
    // Write histogram.
    ofstream fileOut(output_name);
    if(!fileOut){
      cout << "Can't open file " << output_name << endl;
      exit(1);
    }

    for (i = 0; i < hist.GetNumberOfBins(); ++i){
      if (printAll == True){
        fileOut << hist.BinToVal(i) << "," << hist(i) << endl;
      } else if (hist(i) > 0){
        fileOut << hist.BinToVal(i) << "," << hist(i) << endl;
      }
    }
    fileOut.close();
  } else {
    // write to std out.
    for (i = 0; i < hist.GetNumberOfBins(); ++i){
      if (printAll == True){
        cout << hist.BinToVal(i) << "," << hist(i) << endl;
      } else if (hist(i) > 0){
        cout << hist.BinToVal(i) << "," << hist(i) << endl;
      }
    }
  }




  return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

