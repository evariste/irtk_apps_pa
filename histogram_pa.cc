#include <irtkImage.h>
#include <irtkHistogram.h>

char *image_name  = NULL;
char *output_name = NULL;
char *mask_name   = NULL;

void usage()
{
  cerr << "" << endl;
  cerr << " histogram_pa [image] <options> " << endl;
  cerr << " " << endl;
  cerr << " Print the histogram of an image to standard out or to a file." << endl;
  cerr << " By default, empty bins are ignored when printing." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -Tp [value]    : Padding, ignore voxels with intensities below this value." << endl;
  cerr << " -all           : Print contents all bins (not just the non-empty ones)." << endl;
  cerr << " -f [file]      : Write output to a file." << endl;
  cerr << " -mask [file]   : Use a mask to define ROI." << endl;
  cerr << " -bins [count]  : Number of bins for histogram." << endl;
  cerr << " -width [value] : Width of bins." << endl;
  cerr << " " << endl;
  cerr << " N.B. Only one of -bins or -width should be set, not both." << endl;
  cerr << " Setting both values causes the program to quit." << endl;
  cerr << " Setting neither value uses a default of 64 bins and the corresponding" << endl;
  cerr << " bin width for the data." << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkHistogram_1D hist;
  irtkRealImage input;
  irtkRealPixel min = 0, max = 0;
  irtkRealPixel *pPix, *pMask;

  int i, voxels, nBins = 0;
  double padding = -1 * FLT_MAX;
  int ok, printAll = False;
  double width = -1.0;

  // Read arguments
  if (argc < 2){
    usage();
  }

  image_name = argv[1];
  argc--;
  argv++;

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
    if ((ok == False) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      padding = atof(argv[1]);
      argc--;
      argv++;
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
    if ((ok == False) && (strcmp(argv[1], "-width") == 0)){
      argc--;
      argv++;
      width = atof(argv[1]);
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

    if (ok == False){
      cout << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image.
  input.Read(image_name);
  voxels = input.GetNumberOfVoxels();

  // Is there a mask?
  if (mask_name != NULL){
    irtkRealImage mask;
    mask.Read(mask_name);

    pMask = mask.GetPointerToVoxels();
    pPix = input.GetPointerToVoxels();

    for(i = 0; i < voxels; ++i){
      if (*pMask <= 0){
        *pPix = padding;
      }
      ++pPix;
      ++pMask;
    }
  }

  // Initialise histogram.
  input.GetMinMax(&min, &max);

  if (min <= padding){
    // Get min of unpadded voxels
    min = FLT_MAX;
    pPix = input.GetPointerToVoxels();
    for(i = 0; i < voxels; ++i, ++pPix){
      if (*pPix > padding && min > *pPix){
        min = *pPix;
      }
    }
    cerr << "Unpadded ";
  }

  cerr << "min and max : " << min << ", " << max << endl;
  
  hist.PutMin(min);
  hist.PutMax(max);

  if (width > 0){
    if (nBins > 0){
      // Both width and bins have been set.
      cerr << "Both of the parameters -bins and -width have been set ... quitting." << endl;
      exit(1);
    } else {
      // Width only set.
      hist.PutWidth(width);
      cerr << "Histogram assigned width " << width << ", using " << hist.GetNumberOfBins() << " bins." << endl; 
    }
  } else {
    if (nBins > 0){
      // Bins only set.
      hist.PutNumberOfBins(nBins);
      cerr << "Histogram assigned " << nBins << " bins, using width " << hist.GetWidth() << endl;
    } else {
      // Neither set.
      hist.PutNumberOfBins(64);
      cerr << "Histogram uses " << hist.GetNumberOfBins() << " bins with width " << hist.GetWidth() << endl;
    }
  }


  // Populate histogram.
  hist.Reset();
  pPix = input.GetPointerToVoxels();

  for(i = 0; i < voxels; ++i, ++pPix){
    if (*pPix > padding){
      hist.AddSample(*pPix);
    }
  }

  if (output_name != NULL){
    // Write histogram.
    ofstream fileOut(output_name);
    if(!fileOut){
      cerr << "Can't open file " << output_name << endl;
      exit(1);
    }

    fileOut << hist.GetNumberOfBins() << " " << hist.NumberOfSamples() << " " << min << " " << max << endl;
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
    cout << hist.GetNumberOfBins() << " " << hist.NumberOfSamples() << " " << min << " " << max << endl;
    for (i = 0; i < hist.GetNumberOfBins(); ++i){
      if (printAll == True){
        cout << hist.BinToVal(i) << "," << hist(i) << endl;
      } else if (hist(i) > 0){
        cout << hist.BinToVal(i) << "," << hist(i) << endl;
      }
    }
  }

}
