/*
 * concatenate_stacks.cc
 *
 *  Created on: Nov 7, 2013
 *      Author: paulaljabar
 */




#include <irtkImage.h>
#include <irtkFileToImage.h>

char **input_names = NULL;
char *tempStr = NULL;
char *output_name = NULL;


void usage()
{
  cerr << "Usage: " << endl;
  cerr << "concatenate_stacks [output] <options>" << endl;
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "" << endl;
  cerr << "-inputs in1 in2 ... inN : Input stacks" << endl;
  cerr << "-stack_number val       : Use a specific stack's coordinate system, " << endl;
  cerr << "                          where val is chosen from 1 ... N. Default  " << endl;
  cerr << "                          is the central stack (odd N), N/2th stack " << endl;
  cerr << "                          (even N)." << endl;
  cerr << "" << endl;
  cerr << "Assumptions are" << endl;
  cerr << "" << endl;
  cerr << " -  That the orientation of the axes of the input images are all " << endl;
  cerr << "    the same." << endl;
  cerr << "" << endl;
  cerr << " -  That the centres of the image grids all form a straight line " << endl;
  cerr << "    in the patient's coordinate system." << endl;
  cerr << "" << endl;
  cerr << " -  That the line is parallel to the z axis of each of the image " << endl;
  cerr << "    grids" << endl;
  cerr << "" << endl;
  cerr << " -  That we can use a particular stack's axes as the coordinate  " << endl;
  cerr << "    system for the combined data" << endl;
  cerr << "" << endl;
  cerr << " -  That the origin of the new grid coincides with the origin of " << endl;
  cerr << "    the chosen stack." << endl;
  cerr << "" << endl;
  cerr << " -  The axes also coincide and we take as many z-slices as needed " << endl;
  cerr << "    'up' and 'down' to encompass all the data." << endl;
  cerr << "" << endl;
  cerr << " -  That as we loop over the data in memory for each stack, where " << endl;
  cerr << "    one ends, the next one begins." << endl;
  cerr << "" << endl;
  cerr << "NB There may be differences in contrast across different stacks." << endl;
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

      cout << "Input files are: " << endl;
      for (i = 0; i < noOfInputs; i++){
        input_names[i] = argv[1];
        cout << "      " << input_names[i] << endl;
        argc--;
        argv++;
      }
      cout << endl << endl;

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


  irtkImage ** stackImage = new irtkImage*[noOfInputs];


  cout << "Combining following stacks: " << endl;
  for (i = 0; i < noOfInputs; i++)
  {
    cout << input_names[i] << endl;
    stackImage[i] = irtkImage::New(input_names[i]);
  }
  cout << endl;



  irtkImage *targetImg = stackImage[frameStackIndex];

  irtkImageAttributes frameStackAttributes;

  frameStackAttributes = targetImg->GetImageAttributes();



  int slicesBefore, slicesAfter;

  slicesBefore = 0;
  slicesAfter  = 0;
  for (n = 0; n < noOfInputs; n++){

    // Some checking.
    if (stackImage[n]->GetX() != targetImg->GetX() || stackImage[n]->GetY() != targetImg->GetY()){
      cerr << "All stacks must have the same in-plane dimensions." << endl;
      exit(1);
    }

    if (stackImage[n]->GetT() > 1){
      cerr << "Only implemented for 3D volumes. Given volume has at least 4 dimensions: " << input_names[n] << endl;
      exit(1);
    }

    // Main purpose of loop:

    if (n < frameStackIndex){
      slicesBefore += stackImage[n]->GetZ();
    }

    if (n > frameStackIndex){
      slicesAfter += stackImage[n]->GetZ();
    }

    // Collect the origins for the first pair of inputs to see if the direction of the stacks matches
    // the direction of the main z axis. See below.
    if (n == 0){
      stackImage[n]->GetOrigin(x1, y1, z1);
    }

    if (n == 1){
      stackImage[n]->GetOrigin(x2, y2, z2);
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

  switch (targetImg->GetScalarType()) {
  case IRTK_VOXEL_CHAR: {
    targetImg  = new irtkGenericImage<char> (targetAttributes);
  }
  break;
  case IRTK_VOXEL_UNSIGNED_CHAR: {
    targetImg = new irtkGenericImage<unsigned char> (targetAttributes);
  }
  break;
  case IRTK_VOXEL_SHORT: {
    targetImg = new irtkGenericImage<short> (targetAttributes);
  }
  break;
  case IRTK_VOXEL_UNSIGNED_SHORT: {
    targetImg = new irtkGenericImage<unsigned short> (targetAttributes);
  }
  break;
  case IRTK_VOXEL_FLOAT: {
    targetImg = new irtkGenericImage<float> (targetAttributes);
    break;
  }
  case IRTK_VOXEL_DOUBLE: {
    targetImg = new irtkGenericImage<double> (targetAttributes);
    break;
  }
  default:
    cerr << "Unknown voxel type for output format" << endl;
    exit(1);
  }

  targetImg->Print();

  int currZoffset = 0;

  for (n = 0; n < noOfInputs; n++){

    for (k = 0; k < stackImage[n]->GetZ(); ++k){
      for (j = 0; j < stackImage[n]->GetY(); ++j){
        for (i = 0; i < stackImage[n]->GetX(); ++i){
          val = stackImage[n]->GetAsDouble(i,j,k);
          targetImg->PutAsDouble(i, j, k + currZoffset, val);
        }
      }
    }

    currZoffset += stackImage[n]->GetZ();
  }


  targetImg->Write(output_name);


  for (n = 0; n < noOfInputs; n++)
    delete stackImage[n];

  delete targetImg;
  delete [] stackImage;
  delete [] input_names;
}




