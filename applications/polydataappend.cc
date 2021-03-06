#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkAppendPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>

void usage(){
  cerr << "polydataappend [N] [input_1] ... [input_N] [output] <-options>"<<endl;
  cerr << "Options:" << endl;
  cerr << "-lcc         : Keep only largest connected component of each input polydata." << endl;
  cerr << "               Default is to use all components if an input has more than one component" << endl;
  cerr << "-keepscalars : Leave input scalar values untouched." << endl;
  cerr << "               Default is to give uniform scalars to each input that index its position" << endl;
  cout << "               in the argument list." << endl;
  exit(1);
}

int main(int argc, char **argv ){

  if (argc <3 ){
    usage();
  }

  int i, j, noOfPoints, inputCount, regionCount;
  bool ok;
  bool lcc = false;
  bool newScalars = true;

  inputCount = atoi(argv[1]);
  argv++;
  argc--;
  cout << inputCount << " input datasets." << endl;

  char **input_names = new char*[inputCount];
  char *output_name;

  for (i = 0; i < inputCount; ++i){
    input_names[i] = argv[1];
    argv++;
    argc--;
  }

  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-lcc") == 0)){
      argc--;
      argv++;
      lcc = true;
      ok  = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-keepscalars") == 0)){
      argc--;
      argv++;
      newScalars = false;
      ok  = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (lcc == true){
    cout << "Warning! Extracting largest connected components before appending." << endl;
  }

  vtkPolyData** inputPolys;
  vtkPolyDataReader** readers;
  vtkPolyDataConnectivityFilter** connFilter;
  vtkCleanPolyData** cleaners;

  inputPolys = new vtkPolyData*[inputCount];
  readers    = new vtkPolyDataReader*[inputCount];
  connFilter = new vtkPolyDataConnectivityFilter*[inputCount];
  cleaners   = new vtkCleanPolyData*[inputCount];

  vtkPolyDataConnectivityFilter* checkingFilter = vtkPolyDataConnectivityFilter::New();
  checkingFilter->SetExtractionModeToAllRegions();

  vtkFloatArray *scalars = vtkFloatArray::New();

  for (i = 0; i < inputCount; ++i){

    cout << "Reading file " << input_names[i] << endl;

    inputPolys[i] = vtkPolyData::New();
    readers[i]    = vtkPolyDataReader::New();
    connFilter[i] = vtkPolyDataConnectivityFilter::New();
    cleaners[i]   = vtkCleanPolyData::New();

    connFilter[i]->SetExtractionModeToLargestRegion();

    readers[i]->SetFileName(input_names[i]);
    readers[i]->Modified();
    readers[i]->Update();
    inputPolys[i] = readers[i]->GetOutput();

    checkingFilter->SetInputData(inputPolys[i]);
    checkingFilter->Modified();
    checkingFilter->Update();

    noOfPoints = inputPolys[i]->GetNumberOfPoints();

    if (newScalars == true){
      // Add a scalar value for current polydata set.
      scalars->Initialize();
      scalars->SetNumberOfComponents(1);
      scalars->SetNumberOfTuples(noOfPoints);
      scalars->SetName("label");

      for (j = 0; j < noOfPoints; j++){
        scalars->InsertTuple1(j, i);
      }
      inputPolys[i]->GetPointData()->AddArray(scalars);
    }

    if (lcc == true){
      regionCount = checkingFilter->GetNumberOfExtractedRegions();
      if (regionCount > 1){
        cout <<" Warning:: using the largest region of ( ";
        cout << regionCount << " )  of this file" << endl;
      }

      connFilter[i]->SetInputData(inputPolys[i]);
      connFilter[i]->Update();
      cleaners[i]->SetInputData(connFilter[i]->GetOutput());

    } else {
      cleaners[i]->SetInputData(inputPolys[i]);
    }

    cleaners[i]->Modified();
    cleaners[i]->Update();
  }

  // Do the appending.
  vtkAppendPolyData *appendFilter = vtkAppendPolyData::New();
  for (i = 0; i < inputCount; ++i){
    cout << " Adding surface: " << i + 1 <<endl;

    if (i == 0){
      appendFilter->SetInputData(cleaners[i]->GetOutput());
    } else {
      appendFilter->AddInputData(cleaners[i]->GetOutput());
    }
  }

  // Get result.
  vtkPolyData* appended =   vtkPolyData::New();
  appendFilter->SetOutput(appended);
  appendFilter->Update();

  // Recalculate normals
  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SplittingOff();
  normalsFilter->SetInputData(appended);
  normalsFilter->Modified();
  normalsFilter->Update();

  // Write out.
  vtkPolyDataWriter *pd_writer = vtkPolyDataWriter::New();
  cout << "Output file " << output_name << endl;
  pd_writer->SetFileName(output_name);
  pd_writer->SetInputData(normalsFilter->GetOutput());
  pd_writer->Update();
  pd_writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
