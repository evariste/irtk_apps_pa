#if (defined HAS_VTK)

#include <irtkImage.h>
#include <irtkImageFunction.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

char *input_image_name = NULL, *output_name = NULL;
char *input_surface_name = NULL;
char *scalars_name = NULL;

void usage()
{
  cerr << "Usage:  polydataassignscalars [scalarsImage] [surfaceIn] [surfaceOut] <options>" << endl;
  cerr << "" << endl;
  cerr << "Assign scalars to the vertices of surfaceIn.  The scalar assigned" << endl;
  cerr << "is linearly interpolated from [scalarsImage] at each vertex." << endl;
  cerr << "" << endl;
  cerr << "Write the result to surfaceOut." << endl;
  cerr << "Options: " << endl;
  cerr << "-name : Name to give to the assigned scalars." << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, noOfPoints;
  bool ok;

  double xmin, xmax, ymin, ymax, zmin, zmax;
  double pt[3];
  float val;

  if (argc < 4){
    usage();
  }

  // Parse image
  input_image_name  = argv[1];
  argc--;
  argv++;
  input_surface_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-XXX") == 0)){
      argc--;
      argv++;
      // do stuff and maybe argv++ etc.
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
      argc--;
      argv++;
      scalars_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Reading image ..." << endl;
  irtkRealImage *scalarImage = new irtkRealImage(input_image_name);

  vtkPolyData *surface = vtkPolyData::New();

  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_surface_name);
  reader->Update();
  surface = reader->GetOutput();

  noOfPoints= surface->GetNumberOfPoints();

  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(noOfPoints);

  irtkImageFunction *interp = NULL;
  interp = new irtkLinearInterpolateImageFunction;
  interp->SetInput(scalarImage);
  interp->Initialize();

  xmin = -0.5;
  ymin = -0.5;
  zmin = -0.5;
  xmax = scalarImage->GetX()-0.5;
  ymax = scalarImage->GetY()-0.5;
  zmax = scalarImage->GetZ()-0.5;

  cerr << "Assigning scalars ... " << endl;
  for (i = 0; i < noOfPoints; ++i){
    surface->GetPoint(i, pt);
    scalarImage->WorldToImage(pt[0], pt[1], pt[2]);

    //Check that surface does not go outside fov of label image.
    if (pt[0] < xmin || pt[0] > xmax ||
        pt[1] < ymin || pt[1] > ymax ||
        pt[2] < zmin || pt[2] > zmax){
          cerr << "Surface outside bounds of image." << endl;
          exit(1);
    }

    val = interp->Evaluate(pt[0], pt[1], pt[2]);
    scalars->SetTuple1(i,val);
  }

  if (scalars_name == NULL){
    scalars_name = "scalars";
  }

  scalars->Modified();
  scalars->SetName(scalars_name);
  surface->GetPointData()->AddArray(scalars);
  surface->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  return 0;
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

