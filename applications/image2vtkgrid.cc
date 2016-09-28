#ifdef HAS_VTK

///////////////////////////////////////////////////////////////

#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRegularPolygonSource.h>
#include <vtkCubeSource.h>
#include <vtkSphereSource.h>
#include <vtkStructuredPoints.h>

#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridWriter.h>

#include <vtkPolyDataWriter.h>
#include <vtkCleanPolyData.h>
#include <vtkGlyph3D.h>

#include <irtkImage.h>
#include <irtkTransformation.h>

char *input_name = NULL, *output_name = NULL;
char *dof_name = NULL;
void usage(){

  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;
  cout << "" << endl;

  exit(1);
}

int main(int argc, char **argv){

//   int nPts;
  int i, j, k, count;
  bool ok;
  double p[3];
//   double xaxis[3], yaxis[3], zaxis[3];
//   double r;
  int xlast, ylast, zlast;
  int xdim, ydim, zdim;

//   double coord, tmp;

  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  irtkTransformation *transformation;
  if (dof_name != NULL){
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }


  // Read image
  irtkGreyImage image;
  image.Read(input_name);

  xdim = image.GetX();
  ydim = image.GetY();
  zdim = image.GetZ();

  /////////////////////////////////////////

  vtkPoints     *points  = vtkPoints::New();
//   vtkFloatArray *scalars = vtkFloatArray::New();
  vtkPolyData   *output  = vtkPolyData::New();
  output->Allocate();

  xlast = xdim - 1;
  ylast = ydim - 1;
  zlast = zdim - 1;

  int xindices[4] = {0, 1, xlast - 1, xlast};
  int yindices[4] = {0, 1, ylast - 1, ylast};
  int zindices[4] = {0, 1, zlast - 1, zlast};

  count = 0;
  for (i = 0; i < xdim; i++){
    for (j = 0; j < 4; j++){
      for (k = 0; k <4; k++){
        p[0] = i;
        p[1] = yindices[j];
        p[2] = zindices[k];
        image.ImageToWorld(p[0], p[1], p[2]);
        transformation->Transform(p[0], p[1], p[2]);
        points->InsertPoint(count, p);
        ++count;
      }
    }
  }

  for (k = 0; k < zdim; k++){
    for (i = 0; i < 4; i++){
      for (j = 0; j <4; j++){
        p[0] = xindices[i];
        p[1] = yindices[j];
        p[2] = k;
        image.ImageToWorld(p[0], p[1], p[2]);
        transformation->Transform(p[0], p[1], p[2]);
        points->InsertPoint(count, p);
        ++count;
      }
    }
  }

  for (j = 0; j < ydim; j++){
    for (i = 0; i < 4; i++){
      for (k = 0; k < 4; k++){
        p[0] = xindices[i];
        p[1] = j;
        p[2] = zindices[k];
        image.ImageToWorld(p[0], p[1], p[2]);
        transformation->Transform(p[0], p[1], p[2]);
        points->InsertPoint(count, p);
        ++count;
      }
    }
  }

  output->SetPoints(points);
  //  output->GetPointData()->AddArray(scalars);
  output->Modified();


  vtkRegularPolygonSource *polygon = vtkRegularPolygonSource::New();
  polygon->SetNumberOfSides(3);
  polygon->SetRadius(0.3);
  polygon->SetNormal(1, 1, 1);

//   vtkCubeSource *cube = vtkCubeSource::New();
//   r= 0.4;
//   cube->SetBounds(-1*r, r, -1*r, r, -1*r, r);

  //  vtkSphereSource *sphere = vtkSphereSource::New();
  //  sphere->SetRadius(1);

  vtkGlyph3D *glyphs = vtkGlyph3D::New();
  glyphs->SetSourceData(polygon->GetOutput());
  //  glyphs->SetSource(cube->GetOutput());
  //  glyphs->SetSource(sphere->GetOutput());

  glyphs->SetInputData(output);
  glyphs->ScalingOff();

  vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
  cleaner->SetInputData(glyphs->GetOutput());
  cleaner->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInputData(cleaner->GetOutput());
  writer->SetFileName(output_name);
  writer->Write();


//   //////////////////////////////////////////////
//   vtkFloatArray *xcoords = vtkFloatArray::New();
//   vtkFloatArray *ycoords = vtkFloatArray::New();
//   vtkFloatArray *zcoords = vtkFloatArray::New();

//   for (i = 0; i < xdim; ++i){
//     coord = i;
//     tmp = 0;
//     image.ImageToWorld(coord, tmp, tmp);
//     xcoords->InsertNextValue(coord);
//   }

//   for (i = 0; i < ydim; ++i){
//     coord = i;
//     tmp = 0;
//     image.ImageToWorld(tmp, coord, tmp);
//     ycoords->InsertNextValue(coord);
//   }

//   for (i = 0; i < zdim; ++i){
//     coord = i;
//     tmp = 0;
//     image.ImageToWorld(tmp, tmp, coord);
//     zcoords->InsertNextValue(coord);
//   }

//   vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
//   rgrid->SetDimensions(xdim, ydim, zdim);
//   rgrid->SetXCoordinates(xcoords);
//   rgrid->SetYCoordinates(ycoords);
//   rgrid->SetZCoordinates(zcoords);


//   vtkRectilinearGridWriter *rgridWriter = vtkRectilinearGridWriter::New();
//   rgridWriter->SetInputData(rgrid);
//   rgridWriter->SetFileName(output_name);
//   rgridWriter->Write();

//   exit(0);

//   ///////////////////////////////////////////

//   irtkGreyImage dummy(image);

//   vtkStructuredPoints *vtkImage = vtkStructuredPoints::New();
//   xaxis[0] = 1; xaxis[1] = 0; xaxis[2] = 0;
//   yaxis[0] = 0; yaxis[1] = 1; yaxis[2] = 0;
//   zaxis[0] = 0; zaxis[1] = 0; zaxis[2] = 1;
//   dummy.PutPixelSize(1, 1, 1);
//   dummy.PutOrigin(0, 0, 0);
//   dummy.PutOrientation(xaxis, yaxis, zaxis);
//   dummy.ImageToVTK(vtkImage);



}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
