#ifdef HAS_VTK

#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

// Default filenames
char *input_name = NULL;

void usage()
{
  cerr << "Usage: polydatameanradius [polydata] \n" << endl;
  cerr << "Mean radius of points in a poly data set from the centre of\n" << endl;
  cerr << "gravity of all the points." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, nPoints;
  double r, sumR;
  double p[3], cg[3];

  // Check command line
  if (argc < 2){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  nPoints = surface->GetNumberOfPoints();

  if (nPoints < 1){
    cerr << "Polydata in " << input_name << " has no points." << endl;
    exit(1);
  }

  // Find centre of gravity.
  cg[0] = cg[1] = cg[2] = 0;

  for (i = 0; i < nPoints; i++){
    surface->GetPoints()->GetPoint (i, p);
    cg[0] += p[0];
    cg[1] += p[1];
    cg[2] += p[2];
  }

  cg[0] /= nPoints;
  cg[1] /= nPoints;
  cg[2] /= nPoints;

  // Find radii from each point to C of G.
  sumR = 0;
  for (i = 0; i < nPoints; i++){
    surface->GetPoints()->GetPoint (i, p);
    r  = (p[0] - cg[0]) * (p[0] - cg[0]);
    r += (p[1] - cg[1]) * (p[1] - cg[1]);
    r += (p[2] - cg[2]) * (p[2] - cg[2]);
    r = sqrt(r);

    sumR += r;
  }

  cout << "Mean radius : " << sumR / ((double) nPoints) << endl;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

