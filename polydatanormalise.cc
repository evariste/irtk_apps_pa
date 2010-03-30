#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

char *input_name = NULL;
char *output_name = NULL;

void usage() {
  cerr << "" << endl;
  cerr << "Usage: polydatanormalise [input] [output] <options>" << endl; // <options>" << endl;
  cerr << "" << endl;
  cerr
      << "  Centre the polydata and scale its points so that its radius is one."
      << endl;
  cerr << "" << endl;
  cerr << "  Default centre is the origin and the default radius is 1." << endl;
  cerr << "  Centre is defined as the centroid of the points." << endl;
  cerr
      << "  Radius is defined as the maximum distance from the centre over all points."
      << endl;
  cerr << "" << endl;
  cerr << "Options:" << endl;
  cerr << "  -centre x y z   : Alternative centre." << endl;
  cerr << "  -radius r       : Alternative radius." << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv) {
  int i, ok;
  double pt[3];

  double centre[3] = { 0, 0, 0 };

  double radius = 1;

  if (argc < 3) {
    usage();
  }

  // Parse arguments.
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  cerr << "Input        : " << input_name << endl;
  cerr << "Output       : " << output_name << endl;

  while (argc > 1) {
    ok = False;
    if ((!ok) && (strcmp(argv[1], "-centre") == 0)) {
      argc--;
      argv++;
      centre[0] = atof(argv[1]);
      argc--;
      argv++;
      centre[1] = atof(argv[1]);
      argc--;
      argv++;
      centre[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-radius") == 0)) {
      argc--;
      argv++;
      radius = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (!ok) {
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);

  vtkPolyData* input = reader->GetOutput();
  input->Update();

  // the points array
  int noOfPoints = input->GetNumberOfPoints();

  // perform the smoothing
  vtkPoints* pts = input->GetPoints();

  double meanx, meany, meanz;
  double minx, miny, minz;
  double maxx, maxy, maxz;
  double x, y, z;
  double dist, maxDist;

  meanx = meany = meanz = 0.0;
  minx = miny = minz = FLT_MAX;
  maxx = maxy = maxz = -1.0 * FLT_MAX;
  maxDist = -1.0 * FLT_MAX;

  // Loop over surface to get initial info.
  for (i = 0; i < noOfPoints; ++i) {
    pts->GetPoint(i, pt);
    x = pt[0];
    y = pt[1];
    z = pt[2];

    if (minx > x)
      minx = x;
    if (miny > y)
      miny = y;
    if (minz > z)
      minz = z;
    if (maxx < x)
      maxx = x;
    if (maxy < y)
      maxy = y;
    if (maxz < z)
      maxz = z;
    meanx += x;
    meany += y;
    meanz += z;

  }

  meanx /= noOfPoints;
  meany /= noOfPoints;
  meanz /= noOfPoints;

  cout << "Centre    : " << meanx << " " << meany << " " << meanz << endl;
  cout << "min max x :" << minx << " " << maxx << endl;
  cout << "min max y :" << miny << " " << maxy << endl;
  cout << "min max z :" << minz << " " << maxz << endl;

  double meanDistOrig = 0;

  for (i = 0; i < noOfPoints; ++i) {
    pts->GetPoint(i, pt);
    x = pt[0];
    y = pt[1];
    z = pt[2];

    x -= meanx;
    y -= meany;
    z -= meanz;

    meanDistOrig += sqrt(x * x + y * y + z * z);
  }

  meanDistOrig /= noOfPoints;


  if (meanDistOrig <= 0){
    cerr << "Surface seems to have zero or negative radius." << endl;
    exit(1);
  }

  cout << "Original radius : " << meanDistOrig << endl;
  cout << "Original centre : " << meanx << ", " << meany << ", " << meanz << endl;


  // Loop over surface to scale data to new radius.
  for (i = 0; i < noOfPoints; ++i) {
    pts->GetPoint(i, pt);
    x = pt[0];
    y = pt[1];
    z = pt[2];

    x = centre[0] + (x - meanx) * radius / meanDistOrig;
    y = centre[1] + (y - meany) * radius / meanDistOrig;
    z = centre[2] + (z - meanz) * radius / meanDistOrig;

    pts->SetPoint(i, x, y, z);
  }

  input->SetPoints(pts);
  input->Update();

  // save as a vtk file
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(input);
  writer->SetFileName(output_name);
  writer->SetFileTypeToASCII();
  writer->Update();
  writer->Write();

  reader->Delete();
  writer->Delete();
  return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ) {
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

