/*
 * polydatafiles2csv.cc
 *
 *  Created on: Feb 13, 2015
 *      Author: paulaljabar
 */



#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
#include <string.h>
#include <sys/stat.h>

#include <abcdUtils.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

char **inputFilenames  = NULL;
char *outputFileName = NULL;



void usage(std::string name){


  cout << "" << endl;
  cout << "" << endl;
  cout << name <<  " [points.vtk points.vtk .... points.vtk] <options> " << endl;
  cout << "" << endl;
  cout << "Write out the coordinates of the points in the given vtk files in a csv format." << endl;
  cout << "All inputs must have the same number of points." << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "Options:" << endl;
  cout << "" << endl;
  cout << "-file [name]  : Write the output to the specified file. " << endl;
  cout << "                If not specified then written to the terminal's standard output." << endl;
  cout << "" << endl;
  cout << "-verbose      : Say more stuff" << endl;

  exit(1);
}


int main(int argc, char **argv){

  bool ok, verbose;
  double p[3];
  int i, count, j;

  std::string exeName;

  exeName = argv[0];

  if (argc < 2)
    usage(exeName);

  int MAX_FILES = 1000;

  // Parse input file names and values
  char **inputFilenames = new char *[MAX_FILES];

  count = 0;
  while ((argc > 1) && (argv[1][0] != '-')) {
    inputFilenames[count] = argv[1];

    if (! is_vtkPolyDataFile(inputFilenames[count])){
      cerr << "Not a vtk polydata file: " << inputFilenames[count] << endl;
      usage(exeName);
    }

    count++;
    argc--;
    argv++;
  }

  if (count < 1){
    cerr << "No files to read." << endl;
    usage(exeName);
  }

  verbose = false;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-verbose") == 0)){
      argc--;
      argv++;
      verbose = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-file") == 0)){
      argc--;
      argv++;
      outputFileName = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (verbose){
    cerr << "Processing " << count << " files " << endl;
    for (i = 0; i < count; i++)
      cerr << "   " << inputFilenames[i] << endl;
    if (outputFileName != NULL)
      cerr << "Writing output to " << outputFileName << endl;
  }


  vtkSmartPointer<vtkPolyDataReader> reader =
      vtkSmartPointer<vtkPolyDataReader>::New();

  reader->SetFileName(inputFilenames[0]);
  reader->Update();

  vtkSmartPointer<vtkPolyData> inputPoints =
      vtkSmartPointer<vtkPolyData>::New();

  inputPoints = reader->GetOutput();

  int nPts = inputPoints->GetNumberOfPoints();

  FILE *output;

  if (outputFileName == NULL){
    output = stdout;
  } else {
    output = fopen (outputFileName,"w");
  }

  fprintf(output, "file");
  for (j = 0; j < nPts; j++){
    fprintf(output, ",pt-%d-x", j+1);
    fprintf(output, ",pt-%d-y", j+1);
    fprintf(output, ",pt-%d-z", j+1);
  }
  fprintf(output, "\n");

  for (i = 0; i < count; i++){

    reader->SetFileName(inputFilenames[i]);
    reader->Update();

    inputPoints = reader->GetOutput();
    if (inputPoints->GetNumberOfPoints() != nPts){
      cerr << "Input files have differing numbers of points." << endl;
      exit(1);
    }

    fprintf(output, "%s", inputFilenames[i]);
    for (j = 0; j < nPts; j++){
      inputPoints->GetPoint(j, p);
      fprintf(output, ",%f,%f,%f", p[0], p[1], p[2]);
    }
    fprintf(output, "\n");
  }


  if (outputFileName != NULL){
    fclose(output);
  }


}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
