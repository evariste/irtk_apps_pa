
#include <irtkImage.h>
#include <iostream>
#include <fstream>
using namespace std;

char *matrix_name = NULL;
char *csvfile_name = NULL;

void usage()
{
  cerr << " Usage: matrix2csv [matrix] <options>" << endl;
  cerr << " " << endl;
  cerr << " Write the data in a IRTK .mat matrix file to a csv file." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -info          : Just print general information." << endl;
  cerr << " -file file.csv : Output to the csv file specified." << endl;
  cerr << " -transpose     : Transpose rows and columns before writing." << endl;

  exit(1);
}


int main(int argc, char **argv)
{
  int i, j;
  bool ok;
  int infoOnly = false;
  int transpose = false;
  
  if (argc < 3){
    usage();
  }

  // Parse image
  matrix_name = argv[1];
  argc--;
  argv++;
  
  
  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-info") == 0)){
      argc--;
      argv++;
      infoOnly = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-transpose") == 0)){
      argc--;
      argv++;
      transpose = true;
      ok = true;
    }    if ((ok == false) && (strcmp(argv[1], "-file") == 0)){
      argc--;
      argv++;
      csvfile_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }   
      
     if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }
  
  irtkMatrix m;
  
  m.Read(matrix_name);

  if (transpose == true){
    m.Transpose();
  }
  
  if (infoOnly == true){
    cout << "Matrix " << matrix_name << " : " << m.Rows() << "x" << m.Cols();
    if (transpose == true)
      cout << " (after transposition) ";
    cout << endl;
    exit(0);
  }
  
  if (csvfile_name == NULL){
    exit(0);
  }

  ofstream fileout (csvfile_name, ios::out);
  for (i = 0; i < m.Rows(); i++){
    for (j = 0; j < m.Cols(); j++){
      
      fileout << m(i, j);
      if (j < m.Cols() - 1)
        fileout << ",";
    }
    fileout << endl;
  }

  fileout.close();
  
}
