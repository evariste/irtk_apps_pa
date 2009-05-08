#include <irtkImage.h>
#include <irtkMatrix.h>

char *input_name = NULL;
char *output_name = NULL;


void usage(){
  cout << "csv2mat [csvFile] [rows] [cols] [outputFile.mat]" << endl;
  cout << "csv file must be space separated." << endl;

  exit(1);
}

int main(int argc, char **argv){

  if (argc < 5){
    usage();
  }

  int i, j, rows, cols;

  input_name =  argv[1];
  argc--;
  argv++;

  rows = atoi(argv[1]);
  argc--;
  argv++;

  cols = atoi(argv[1]);
  argc--;
  argv++;

  output_name =  argv[1];
  argc--;
  argv++;

  irtkMatrix m;
  m.Initialize(rows, cols);

  ifstream in(input_name);

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      in >> m(i, j);
    }
  }

  m.Print();

  m.Write(output_name);

}
