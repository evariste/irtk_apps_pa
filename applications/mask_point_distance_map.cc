#include <irtkImage.h>

#include <abcdPointDistanceMap.h>

char *input_name = NULL, *output_name = NULL;
char *pointFileName = NULL;

void usage()
{
  cerr << "Usage: mask_point_distance_map [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t -pointsFile file   VTK format file containing world coordinate location(s)" << endl;
  cerr << "\t -point_i a b c     Image coordinates of a single point."<< endl;
  cerr << "\t -point_w x y z     World coordinates of a single point."<< endl;

  exit(1);
}


int main(int argc, char *argv[])
{

  int ok, i;
  double pt[3];

  // Check command line
  if (argc < 3) {
    usage();
  }



  // Read input and output names
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;


  abcdPointDistanceMap<irtkGreyPixel> ptDmap;


  irtkGreyImage *input = new irtkGreyImage(input_name);
  irtkGreyImage *output = new irtkGreyImage(input_name);

  ptDmap.SetInput(input);
  ptDmap.SetOutput(output);


  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-pointsFile") == 0)) {
      argc--;
      argv++;
      pointFileName = argv[1];
      argc--;
      argv++;
      ok = true;

    }
    if ((ok == false) && (strcmp(argv[1], "-point_i") == 0)) {
      argc--;
      argv++;
      ptDmap.addSeedPointI(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));

      argc -= 3;
      argv += 3;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-point_w") == 0)) {
      argc--;
      argv++;
      ptDmap.addSeedPointW(atof(argv[1]), atof(argv[2]), atof(argv[3]));

      argc -= 3;
      argv += 3;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  if (pointFileName != NULL){
    vtkPolyDataReader *pd_reader = vtkPolyDataReader::New();
    pd_reader->SetFileName(pointFileName);
    pd_reader->Modified();
    pd_reader->Update();
    vtkPolyData *pd = pd_reader->GetOutput();
    for (i = 0; i < pd->GetNumberOfPoints(); ++i){
      pd->GetPoint(i, pt);
      ptDmap.addSeedPointW(pt[0], pt[1], pt[2]);
    }
  }

  if (ptDmap.GetNumberOfSeedPoints() < 1){
    cerr << "mask_point_distance_map : Need at least one seed point" << endl;
    exit(1);
  }


  ptDmap.Run();

	output->Write(output_name);


}

