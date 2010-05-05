#include <irtkImage.h>

#include <irtkHistogram.h>

#include <irtkTransformation.h>

// Default filenames
char *output_name = NULL, **input_name = NULL;
char *sd_name = NULL;

void usage()
{
  cerr << "Usage: atlas-jacobians [output] [N] [input1..inputN] <options>\n" 
       << endl;
  cerr << "Calculate voxel-wise geometric mean (via logs) of jacobian images (scaled as a percentage growth)." << endl;
  cerr << "where <options> is one or more of the following:\n";
  cerr << "<-scale factor>  Scaling factor (0 means no normalsation)\n";
  cerr << "<-sd file>       Write standard deviation of logs to file." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double scale, val;
  int i, j, n, padding, no, ok;

  irtkRealPixel *ptr2input;
  irtkRealImage input;

  irtkRealPixel *ptr2mean;
  irtkRealImage meanLog;

  irtkRealPixel *ptr2sd;
  irtkRealImage sdLog;

  irtkGreyImage counts;
  irtkGreyPixel *ptr2counts;

  // Check command line
  if (argc < 5){
    usage();
  }

  // Parse image
  output_name = argv[1];
  argc--;
  argv++;
  
  // Parse number of input images
  no = atoi(argv[1]);
  argc--;
  argv++;

  // Parse input file names
  input_name = new char *[no];
  for (i = 0; i < no; i++){
    input_name[i] = argv[1];
    argc--;
    argv++;
  }

  // Default scaling factor
  scale = 0;

  // Default padding 
  padding = MIN_GREY;

  // Parse any remaining paramters
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-scale") == 0)){
      argc--;
      argv++;
      scale = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sd") == 0)){
      argc--;
      argv++;
      sd_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 
  
  // Read and add one input image after the other
  n = 0;
  for (i = 0; i < no; i++){
    
    // Read input
    cout << "Reading input image " << input_name[i] << "... "; cout.flush();
    input.Read(input_name[i]);
    cout << "done" << endl;

    if (i == 0){
      n = input.GetX()*input.GetY()*input.GetZ();
      meanLog = input;
      sdLog   = input;
      counts  = input;

      ptr2mean   = meanLog.GetPointerToVoxels();
      ptr2sd     = sdLog.GetPointerToVoxels();
      ptr2counts = counts.GetPointerToVoxels();

      for (j = 0; j < n; j++){
	*ptr2mean = 0;
        *ptr2sd   = 0;
        *ptr2counts = 0;

	++ptr2mean;
        ++ptr2sd;
        ++ptr2counts;
      }
    }

    if (scale > 0){

      cout << "Not implemented with scale ... " << endl;
      exit(1);
      
    } else {

      cerr << "Adding input to atlas..."; cout.flush();

      ptr2input  = input.GetPointerToVoxels();
      ptr2mean   = meanLog.GetPointerToVoxels();
      ptr2sd     = sdLog.GetPointerToVoxels();
      ptr2counts = counts.GetPointerToVoxels();

      for (j = 0; j < n; j++){

        if (*ptr2input > 0){
          val = log(*ptr2input / 100.0);
          *ptr2mean += val;
          *ptr2sd   += val * val;
          *ptr2counts += 1;
        }

	++ptr2input;
	++ptr2mean;
        ++ptr2sd;
        ++ptr2counts;
      }
      cerr << "done\n";

    }
  } 

  cerr << "Scaling atlas by factor " << scale << " ... "; 
  cout.flush();
  if (scale > 0){

  } else {

    ptr2input  = input.GetPointerToVoxels();
    ptr2mean   = meanLog.GetPointerToVoxels();
    ptr2sd     = sdLog.GetPointerToVoxels();
    ptr2counts = counts.GetPointerToVoxels();

    for (j = 0; j < n; j++){

      if (*ptr2counts > 0){
        val = *ptr2mean / (double)(*ptr2counts);
        *ptr2input = 100.0 * exp (val) ;
        *ptr2sd = sqrt((*ptr2sd / (double)(*ptr2counts)) - val * val);
      } else {
        *ptr2input = 0.0;
      }

      ++ptr2input;
      ++ptr2mean;
      ++ptr2sd;
      ptr2counts++;
    }
  }

  cerr << " done\n";

  // Write atlas
  cerr << "Writing atlas to " << output_name << " ... "; 
  cout.flush();
  input.Write(output_name);
  cerr << " done" << endl;

  if (sd_name != NULL){
    cerr << "Writing SD of logs to " << sd_name << " ...";
    cout.flush();
    sdLog.Write(sd_name);
    cerr << " done" << endl;
  }

}
