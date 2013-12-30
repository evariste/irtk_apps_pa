

#include <irtkImage.h>
#include <irtkFileToImage.h>

char **input_names = NULL;
char *tempStr = NULL;
char *output_name = NULL;


void usage()
{
  cerr << "Usage: " << endl;
  cerr << "NAME_OF_EXE_HERE [output] <options>" << endl;
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "" << endl;
  cerr << "-inputs input1 input2 ... inputN      Input stacks" << endl;
  cerr << "-stack_number val                     Use a specific stack's coordinate system, " << endl;
  cerr << "                                      where val is chosen from 1 ... N. Default is " << endl;
  cerr << "                                      the central stack (odd N), N/2th stack (even N)." << endl;
  cerr << "" << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  int ok, i, j, k, n, noOfInputs, frameStackIndex;
  bool readingInput;
  double x1, y1, z1, x2, y2, z2, val;


  noOfInputs = 0;
  frameStackIndex = -1;


  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-inputs") == 0)) {
      argc--;
      argv++;

      readingInput = true;
      while (noOfInputs < argc - 1 && readingInput){
        noOfInputs += 1;
        if (argv[noOfInputs][0] == '-'){
          readingInput = false;
          noOfInputs--;
        }
      }

      input_names = new char*[noOfInputs];

      for (i = 0; i < noOfInputs; i++){
        input_names[i] = argv[1];
        argc--;
        argv++;
      }
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-stack_number") == 0)){
      argc--;
      argv++;

      frameStackIndex = atoi(argv[1]) - 1;

      argc--;
      argv++;

      ok = true;
    }


    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  if (noOfInputs < 2){
    cerr << "At least 2 inputs should be given." << endl;
    exit(1);
  }

  if (frameStackIndex < 0){
    frameStackIndex = noOfInputs / 2;
  }

  if (frameStackIndex > noOfInputs - 1){
    cerr << "Stack number out of range: " << frameStackIndex + 1 << endl;
    cerr << "Exceeds number of inputs." << endl;
    frameStackIndex = noOfInputs / 2;
  }

  cout << "Using frame from stack number " << frameStackIndex + 1 << ": " << input_names[frameStackIndex] << endl;


  cout << "Combining following stacks: " << endl;
  for (i = 0; i < noOfInputs; i++)
    cout << input_names[i] << endl;
  cout << endl;

  irtkFileToImage *frameStackReader = irtkFileToImage::New(input_names[frameStackIndex]);
  irtkBaseImage *targetImg;
  targetImg = frameStackReader->GetOutput();

  irtkImageAttributes frameStackAttributes;

  frameStackAttributes = targetImg->GetImageAttributes();

  /*
   * Assumptions are that the orientation of the axes of the input images are all the same.
   * That the centres of the image grids all form a straight line in the patients coordinate system.
   * That the line is parallel
   * to the z axis of each of the image grids.
   * That we can use the one of the image grids as a coordinate system for the
   * combined data
   * That the origin of the new grid coincides with the origin of the chosen stack.
   * The axes also coincide and we take as many z-slices as needed 'up' and 'down' to encompass all the data.
   * That as we loop over the data in memory for each stack, where one ends, the next one begins.
   *
   * There may be differences in contrast across different stacks.
   */

  int slicesBefore, slicesAfter;

  slicesBefore = 0;
  slicesAfter  = 0;
  for (n = 0; n < noOfInputs; n++){

    irtkFileToImage *reader = irtkFileToImage::New(input_names[n]);
    irtkBaseImage *image = reader->GetOutput();

    // Some checking.
    if (image->GetX() != targetImg->GetX() || image->GetY() != targetImg->GetY()){
      cerr << "All stacks must have the same in-plane dimensions." << endl;
      exit(1);
    }

    if (image->GetT() > 1){
      cerr << "Only implemented for 3D volumes. Given volume has at least 4 dimensions: " << input_names[n] << endl;
      exit(1);
    }

    // Main purpose of loop:

    if (n < frameStackIndex){
      slicesBefore += image->GetZ();
    }

    if (n > frameStackIndex){
      slicesAfter += image->GetZ();
    }

    // Collect the origins for the first pair of inputs to see if the direction of the stacks matches
    // the direction of the main z axis. See below.
    if (n == 0){
      image->GetOrigin(x1, y1, z1);
    }

    if (n == 1){
      image->GetOrigin(x2, y2, z2);
    }

  }

  // Direction of stacks as estimated by the displacement of the origins from first to last:
  x1 = x2 - x1;
  y1 = y2 - y1;
  z1 = z2 - z1;

  // Inner product of the estimated direction with the z-axis of the stack used as a frame.
  val = x1 * frameStackAttributes._zaxis[0] + y1 * frameStackAttributes._zaxis[1] + z1 * frameStackAttributes._zaxis[2];

  if (val < 0){
    // Stacks are ordered in opposite direction to z axis.
    swap(slicesAfter, slicesBefore);
  }

  // Offset of origin.

  double dispX, dispY, dispZ;

  // z-axis components:
  x1 = frameStackAttributes._zaxis[0];
  y1 = frameStackAttributes._zaxis[1];
  z1 = frameStackAttributes._zaxis[2];

  val = sqrt(x1*x1 + y1*y1 + z1*z1);

  if (val <= 0){
    cerr << "Value error in norm of z axis." << endl;
    exit(1);
  }

  x1 /= val;
  y1 /= val;
  z1 /= val;


  // Required offset in mm.
  val = 0.5 * (slicesAfter - slicesBefore) * targetImg->GetZSize();

  dispX = val * x1;
  dispY = val * y1;
  dispZ = val * z1;


  frameStackAttributes._xorigin += dispX;
  frameStackAttributes._yorigin += dispY;
  frameStackAttributes._zorigin += dispZ;

  int zdimTotal;

  zdimTotal = slicesBefore + slicesAfter + targetImg->GetZ();

  irtkImageAttributes targetAttributes = frameStackAttributes;

  targetAttributes._z = zdimTotal;
  targetImg->Initialize(targetAttributes);
  targetImg->Print();

  int currZoffset = 0;

  for (n = 0; n < noOfInputs; n++){
    irtkFileToImage *reader = irtkFileToImage::New(input_names[n]);
    irtkBaseImage *image = reader->GetOutput();

    for (k = 0; k < image->GetZ(); ++k){
      for (j = 0; j < image->GetY(); ++j){
        for (i = 0; i < image->GetX(); ++i){
          val = image->GetAsDouble(i,j,k);
          targetImg->PutAsDouble(i, j, k + currZoffset, val);
        }
      }
    }

    currZoffset += image->GetZ();
  }



  targetImg->Write(output_name);

}
