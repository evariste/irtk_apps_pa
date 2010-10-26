/*=========================================================================

  Library   : packages/segmentation
  Module    : $RCSfile: ems.cc,v $
  Authors   : Maria Murgasova and Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2006
  Purpose   :
  Date      : $Date: 2008-06-05 09:38:41 $
  Version   : $Revision: 1.7 $
  Changes   : $Locker:  $
              $Log: ems.cc,v $
              Revision 1.7  2008-06-05 09:38:41  mm3
              added option for not adding extra tissue for background

              Revision 1.6  2008-05-19 16:11:31  mm3
              changed to irtkRealImage, added app bci

              Revision 1.5  2007-05-24 14:48:22  mm3
              saving probability maps dor 7 structure
              Save probability maps for seven structure.

              Revision 1.4  2007/02/15 10:37:23  dr
              Fixed missing bracket

              Revision 1.3  2007/02/12 12:12:33  mm3
              Write probmaps for 3 or 11 tissues

              Revision 1.2  2007/02/12 11:52:21  mm3
              option -background added
              writes the probability maps into current directory

              Revision 1.1  2006/10/15 21:15:10  dr
              Imported sources


=========================================================================*/

#include <irtkImage.h>

#include <irtkEMClassification.h>

#include <irtkGaussian.h>

char *output_name;

char *soft_csf_name   = NULL;
char *soft_wm_name    = NULL;
char *soft_gm_name    = NULL;
char *parameters_name = NULL;
char *mask_name       = NULL;

void usage()
{
  cout << "Usage: ems [image] [n] [atlas 1 ... atlas n] [output] <options>" << endl;
  cout << " Where options are : " << endl;
  cout << " -iterations [val]" << endl;
  cout << " -padding [val]            : Padding value." << endl;
  cout << " -mask [file]              : Mask that to use to define ROI (>0)." << endl;
  cout << " -background [file]        : Background probability map." << endl;
  cout << " -softmaps [csf] [gm] [wm] : Output soft segmentations to these files." << endl;
  cout << " -parameter file           : File to write parameters to." << endl;
  cout << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, padding, iterations;
  bool ok, nobg=false;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  image.Read(argv[1]);
  argc--;
  argv++;

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;
  
  // Probabilistic atlas
  irtkRealImage **atlas = new irtkRealImage*[n];
  irtkRealImage *background=NULL;

  // Read atlas for each tissue
  for (i = 0; i < n; i++){
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    cout << "Image " << i <<"            = " << argv[1] <<endl;
    argc--;
    argv++;
  }

  // File name for output
  output_name = argv[1];
  argc--;
  argv++;

  // Default parameters
  iterations = 15;
  padding    = -1;//MIN_GREY;

  // Parse remaining parameters
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
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-background") == 0)){
      argc--;
      argv++;
      background = new irtkRealImage;
      background->Read(argv[1]);
      cout << "Background map     = " << argv[1] <<endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-parameter") == 0)){
      argc--;
      argv++;
      parameters_name = argv[1];
      cout << "Parameters file    = " << parameters_name <<endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
      argc--;
      argv++;
      mask_name = argv[1];
      cout << "Using mask            " << mask_name <<endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-softmaps") == 0)){
      argc--;
      argv++;
      soft_csf_name = argv[1];
      argc--;
      argv++;
      soft_gm_name = argv[1];
      argc--;
      argv++;
      soft_wm_name = argv[1];
      argc--;
      argv++;
      ok = true;
      cout << "Writing soft to csf    " << soft_csf_name << endl;
      cout << "Writing soft to wm     " << soft_wm_name << endl;
      cout << "Writing soft to gm     " << soft_gm_name << endl;
    }
    if (ok == false){
      cout << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  if (mask_name != NULL){
    irtkGreyImage mask;
    mask.Read(mask_name);

    int voxels = image.GetNumberOfVoxels();
    if (voxels != mask.GetNumberOfVoxels()){
      cerr << "Mask and image have different sizes." << endl;
      exit(1);
    }
    cout << "Applying mask ... ";
    irtkRealPixel *ptr2image = image.GetPointerToVoxels();
    irtkGreyPixel *ptr2mask  = mask.GetPointerToVoxels();
    for (i =0; i < voxels; ++i){
      if (*ptr2mask < 1){
        *ptr2image = padding;
      }
      ++ptr2mask;
      ++ptr2image;
    }
    cout << "done." << endl;
  }


  irtkEMClassification classification(n, atlas, background);

  classification.SetInput(image);
  classification.SetPadding(padding);
  classification.Initialise();

  double rel_diff;
  //for (i = 0; rel_diff > 0.0001; i++){
  i=0;
  cout << "------------------------------------" << endl;
  do
  {
    cout << "Iteration = " << i+1 << " / " << iterations << endl;
    rel_diff = classification.Iterate(i);
    i++;
    cout << "------------------------------------" << endl;
  }
  while((rel_diff > 0.001) && (i < 50) && (i < iterations));


  irtkRealImage segmentation;
  classification.ConstructSegmentation(segmentation);
  segmentation.Write(output_name);
  if(n==11)
  {
    classification.WriteProbMap(0,"csf.hdr");
    classification.WriteProbMap(1,"gray.hdr");
    classification.WriteProbMap(2,"caudate.hdr");
    classification.WriteProbMap(3,"putamen.hdr");
    classification.WriteProbMap(4,"nigra.hdr");
    classification.WriteProbMap(5,"cerebellum.hdr");
    classification.WriteProbMap(6,"thalamus.hdr");
    classification.WriteProbMap(7,"pallidum.hdr");
    classification.WriteProbMap(8,"brainstem.hdr");
    classification.WriteProbMap(9,"white.hdr");
    classification.WriteProbMap(10,"cerebellum-white.hdr");
    classification.WriteProbMap(11,"other.hdr");
  }
  if(n==7)
  {
    classification.WriteProbMap(0,"csf.hdr");
    classification.WriteProbMap(1,"gray.hdr");
    classification.WriteProbMap(2,"caudate.hdr");
    classification.WriteProbMap(3,"putamen.hdr");
    classification.WriteProbMap(4,"thalamus.hdr");
    classification.WriteProbMap(5,"pallidum.hdr");
    classification.WriteProbMap(6,"white.hdr");
    classification.WriteProbMap(7,"other.hdr");
  }
  if (n==3 && soft_csf_name != NULL && soft_gm_name != NULL && soft_wm_name != NULL)
  {
    classification.WriteProbMap(0, soft_csf_name);
    classification.WriteProbMap(1, soft_gm_name);
    classification.WriteProbMap(2, soft_wm_name);
  }

  if (parameters_name == NULL){
    classification.WriteGaussianParameters("parameters.txt");
  } else {
    classification.WriteGaussianParameters(parameters_name);
  }

  
}

