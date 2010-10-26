#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "  Usage: addslices [in] [out] <-x xslices> <-y yslices> <-z zslices>" << endl;
  cerr << "  Add slices symmetrically on either side of an image. E.g.:" << endl;
  cerr << "  " << endl;
  cerr << "   addslices in.nii.gz out.nii.gz -x 3" << endl;
  cerr << "" << endl;
  cerr << "  adds 6 slices in the x direction 3 on either side of the original volume." << endl;
  cerr << "  The new slices have the same voxel dimensions as the original (of course)." << endl;
  cerr << "  " << endl;
  cerr << "  Options:  -real    Write output as floating point image." << endl;
  cerr << "  " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  irtkRealImage in;

  int xslices = 0, yslices = 0, zslices = 0;
  int xdim, ydim, zdim;
  int xdimnew, ydimnew, zdimnew;
  int i, j, k;

  bool useReal = false;

  double xsize, ysize, zsize;
  double oX, oY, oZ;
  irtkPoint origin;
  double xaxis[3];
  double yaxis[3];
  double zaxis[3];

  if (argc < 3){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      xslices = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      yslices = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-z") == 0)){
      argc--;
      argv++;
      zslices = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-real") == 0)){
      argc--;
      argv++;
      useReal = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }


  // Read input
  in.Read(input_name);

//   if (argc == 0){
//     cout << "Nothing to add, exiting." << endl;
//     exit(0);
//   }




  xdim = in.GetX();
  ydim = in.GetY();
  zdim = in.GetZ();

  xdimnew = 2 * xslices + xdim;
  ydimnew = 2 * yslices + ydim;
  zdimnew = 2 * zslices + zdim;

  xsize = in.GetXSize();
  ysize = in.GetYSize();
  zsize = in.GetZSize();

  in.GetOrigin(oX, oY, oZ);
  origin._x = oX;
  origin._y = oY;
  origin._z = oZ;

  in.GetOrientation(&xaxis[0], &yaxis[0], &zaxis[0]);

  irtkImageAttributes attr;
  attr._x = xdimnew;
  attr._y = ydimnew;
  attr._z = zdimnew;

  attr._dx = xsize;
  attr._dy = ysize;
  attr._dz = zsize;

  attr._xorigin = origin._x;
  attr._yorigin = origin._y;
  attr._zorigin = origin._z;

  for (i = 0; i < 3; ++i){
    attr._xaxis[i] = xaxis[i];
    attr._yaxis[i] = yaxis[i];
    attr._zaxis[i] = zaxis[i];
  }

  irtkRealImage out(attr);

//  irtkRealImage out(xdimnew, ydimnew, zdimnew, xsize, ysize, zsize,
//                   origin, xaxis, yaxis, zaxis);

//   out->PutRegion(xslices, yslices, zslices,
//                 xdim + xslices - 1,
//                 zdim + yslices - 1,
//                 zdim + zslices - 1,
//                 in);

  cout << xslices << " " << yslices << " " << zslices << " " << endl;

  out.Print();


  for (k = 0; k < zdim; ++k){
    for (j = 0; j < ydim; ++j){
      for (i = 0; i < xdim; ++i){
        out.Put(i + xslices, j + yslices, k + zslices, in.Get(i, j, k));
      }
    }
  }


  if (useReal == true){
    cerr << "Writing floating point image." << endl;
    out.Write(output_name);
  } else {
    irtkGreyImage out2(out);
    out2.Write(output_name);
  }

}

