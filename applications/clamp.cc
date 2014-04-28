#include <irtkImage.h>
#include <nr.h>


char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "\tclamp [input] [output] <options>" << endl;
  cerr << "\tClamp the values in [input] according to percentages specified." << endl;
  cerr << "\tE.g if 2 (%) is required for both upper and lower ranges, the values" << endl;
  cerr << "\tare clamped to the 2nd and 98th percentiles in the image." << endl;
  cerr << "\tOptions:" << endl;
  cerr << "\t-p val : percentage to be applied to both upper and lower ends of range." << endl;
  cerr << "\t         In this case lower is val and upper is 100 - val." << endl;
  cerr << "\t-l val : lower percentage value. (default 0)" << endl;
  cerr << "\t-u val : upper percentage value. (default 100)" << endl;
  cerr << "\t-pad value : ignore voxels <= value." << endl;
  cerr << "\t" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 4){
    usage();
  }

  irtkRealPixel *ptr;

  double pad = -1.0 * FLT_MAX;
  double lowerPercentage = 0;
  double upperPercentage = 100;

  int count, voxels, i;
  bool ok;
  float *values;
  int minIndex, maxIndex;
  double minVal, maxVal;

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;      argv++;
      pad = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-p") == 0)){
      argc--;      argv++;
      lowerPercentage = atof(argv[1]);
      upperPercentage = 100.0 - atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-l") == 0)){
      argc--;      argv++;
      lowerPercentage = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-u") == 0)){
      argc--;      argv++;
      upperPercentage = atof(argv[1]);
      argc--;      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  irtkRealImage *input = new irtkRealImage(input_name);

  // How many unpadded voxels are there?
  voxels = input->GetNumberOfVoxels();
  ptr = input->GetPointerToVoxels();
  count = 0;

  for (i = 0; i < voxels; ++i){
    if (*ptr > pad){
      ++count;
    }
    ++ptr;
  }

  values = new float[count + 1];


  ptr = input->GetPointerToVoxels();
  count = 0;
  for (i = 0; i < voxels; ++i){
    if (*ptr > pad){
      ++count;
      values[count] = *ptr;
    }
    ++ptr;
  }

  sort(count, values);

  minIndex = (int) round (count * lowerPercentage / 100.0);
  maxIndex = (int) round (count * upperPercentage / 100.0);

  minVal = values[minIndex];
  maxVal = values[maxIndex];

  cout << "Using min and max : " << minVal << " " << maxVal << " for clamping" << endl;

  ptr = input->GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    if (*ptr > pad){
      if (*ptr > maxVal)
        *ptr = maxVal;

      if (*ptr < minVal)
        *ptr = minVal;
    }
    ++ptr;
  }

  input->Write(output_name);

  delete [] values;

}



