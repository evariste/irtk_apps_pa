#include <irtkImage.h>

//#include <nr.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>


// Default filenames
char *output_name = NULL, **input_name = NULL;

#define MAX_IMAGES 100

void usage()
{
  cerr << "Usage: imagesmedian [output] [input1..inputN]\n"
       << endl;
  cerr << "";
  cerr << "";
  exit(1);
}


int main(int argc, char **argv)
{
  int i, no, j, k, n;
  int zdim, ydim, xdim;

  if (argc < 3)
    usage();

  output_name = argv[1];
  argc--;
  argv++;

  input_name = new char *[MAX_IMAGES];

  // Parse any remaining paramters
  i = 0;
  while (argc > 1){
    input_name[i] = argv[1];
    argc--;
    argv++;
    i++;
  }
  no = i;

  irtkGreyImage *input;
  irtkGreyImage output;

  output.Read(input_name[0]);

  input = new irtkGreyImage[no];
  for (i = 0; i < no; ++i){
    input[i].Read(input_name[i]);
  }

//  float *data;
//  data = new float[1 + no];
//  data[0] = 0;
  gsl_vector *data = gsl_vector_alloc(no);


  xdim = input[0].GetX();
  ydim = input[0].GetY();
  zdim = input[0].GetZ();

  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){

        for (n = 0; n < no; ++n){
//          data[1 + n] = input[n].Get(i, j, k);
          gsl_vector_set(data, n, input[n].Get(i, j, k));

        }

//        sort(no, data);
        gsl_sort_vector(data);

//        output.Put(i, j, k, (int) data[(1 + no)/2]);
        output.Put(i, j, k, (int) gsl_vector_get(data, no/2));

      }
    }
  }

  output.Write(output_name);


}
