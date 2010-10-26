#include <irtkImage.h>

#include <irtkDilationWithPadding.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: dilationWithPadding [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-iterations n>    Number of iterations\n";
  cerr << "\t<-pad value>       Padding\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int i, iterations;
  irtkGreyImage image;
  irtkGreyPixel pad = 0;
 
  // Check command line
  if (argc < 3){
    usage();
  }

  // Read input and output names
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  image.Read(input_name);

  // Parse remaining parameters
  iterations = 1;
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)){
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;
      argv++;
      pad = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  cout << "Dilating with padding value " << pad << " ... "; cout.flush();
  irtkDilationWithPadding<irtkGreyPixel> dilation;
  dilation.PutPaddingValue(pad);
  dilation.SetInput(&image);
  dilation.SetOutput(&image);
  for (i = 0; i < iterations; i++) dilation.Run();
  cout << "done" << endl;

  // Save result
  image.Write(output_name);

  return 0;
}


