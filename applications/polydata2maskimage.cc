#ifdef HAS_VTK


#include <irtkImage.h>

#include <irtkGaussianBlurring.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkStripper.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencil.h>
#include <vtkImageData.h>
#include <vtkTriangleFilter.h>

char *input_polys = NULL, *output_name = NULL;
char *template_image_name = NULL;

void usage()
{
  cerr << "polydata2maskimage [polydata] [template image] [output image] <options>" << endl;
  cerr << "-value val : value to put into the structure (default 1000)." << endl;
  cerr << "-reverse   : toggle whether the structure or the background gets the value." << endl;
  cerr << "" << endl;
  cerr << "See also: http://www.cmake.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataToImageData" << endl;
  cerr << "for similar code." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int n, i, j, k, xdim, ydim, zdim;
  irtkGreyPixel *ptr1;
  short *ptr2;
  short val;
  int noOfPts;
  double pt[3], newPt[3];
  double bounds[6];
  int i1, j1, k1, i2, j2, k2;
  double value = 1000;
  int reverse = true;
  double xspacing, yspacing, zspacing;


  if (argc < 4){
    usage();
  }

  // Parse parameters
  input_polys = argv[1];
  argv++;
  argc--;
  template_image_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-value") == 0)){
      argc--;
      argv++;
      value = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-reverse") == 0)){
      argc--;
      argv++;
      reverse = true;
      ok = true;
    }

    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_polys);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *inputRead = surface_reader->GetOutput();

  vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
  triFilter->SetInput(inputRead);
  triFilter->Update();


  vtkPolyData* polys = vtkPolyData::New();
  polys = triFilter->GetOutput();
  polys->Update();



  // Read image
  irtkGreyImage image, imageSubregion;
  image.Read(template_image_name);

  noOfPts = polys->GetNumberOfPoints();

  // What are the extents of the polydata in image coordinates?
  i1 = image.GetX() - 1;
  i2 = 0;
  j1 = image.GetY() - 1;
  j2 = 0;
  k1 = image.GetZ() - 1;
  k2 = 0;

  for (i = 0; i < noOfPts; ++i){
    polys->GetPoints()->GetPoint (i, pt);
    image.WorldToImage(pt[0], pt[1], pt[2]);
    if (i1 > pt[0])
      i1 = (int) floor(pt[0]) - 1;
    if (j1 > pt[1])
      j1 = (int) floor(pt[1]) - 1;
    if (k1 > pt[2])
      k1 = (int) floor(pt[2]) - 1;
    if (i2 < pt[0])
      i2 = (int) floor(pt[0]) + 2;
    if (j2 < pt[1])
      j2 = (int) floor(pt[1]) + 2;
    if (k2 < pt[2])
      k2 = (int) floor(pt[2]) + 2;
  }

  // Clamp.
  i1 = max(i1, 0);
  j1 = max(j1, 0);
  k1 = max(k1, 0);
  i2 = min(i2, image.GetX() - 1);
  j2 = min(j2, image.GetY() - 1);
  k2 = min(k2, image.GetZ() - 1);

  imageSubregion = image.GetRegion(i1, j1, k1, i2, j2, k2);

  irtkMatrix w2i = imageSubregion.GetWorldToImageMatrix();

  xdim = imageSubregion.GetX();
  ydim = imageSubregion.GetY();
  zdim = imageSubregion.GetZ();

  // Create a blank image where world and image coordinates coincide.
  irtkGreyImage imageCanonical(xdim, ydim, zdim);
  imageCanonical.PutOrigin((xdim - 1)/2.0, (ydim - 1)/2.0, (zdim - 1)/2.0);

  ptr1 = imageCanonical.GetPointerToVoxels();
  n = imageCanonical.GetNumberOfVoxels();
  for (i = 0; i < n; ++i){
    *ptr1 = 0;
    ++ptr1;
  }

  // Convert the polydata points to canonical image coordinates.

  for (i = 0; i < noOfPts; ++i){
    polys->GetPoints()->GetPoint (i, pt);
    newPt[0] = w2i(0, 0) * pt[0] + w2i(0, 1) * pt[1] + w2i(0, 2) * pt[2] + w2i(0, 3);
    newPt[1] = w2i(1, 0) * pt[0] + w2i(1, 1) * pt[1] + w2i(1, 2) * pt[2] + w2i(1, 3);
    newPt[2] = w2i(2, 0) * pt[0] + w2i(2, 1) * pt[1] + w2i(2, 2) * pt[2] + w2i(2, 3);
    polys->GetPoints()->SetPoint(i, newPt);
  }
  polys->Modified();
  polys->Update();

  // Convert image to VTK format to give to stencil.
  vtkStructuredPoints *vtkimageIn = vtkStructuredPoints::New();
  imageCanonical.ImageToVTK(vtkimageIn);

  vtkPolyDataToImageStencil *dataToStencil = vtkPolyDataToImageStencil::New();
  dataToStencil->SetTolerance(0.0);
  dataToStencil->SetInput(polys);

  dataToStencil->SetInformationInput(vtkimageIn);


  vtkImageStencil *stencil = vtkImageStencil::New();
  stencil->SetInput(vtkimageIn);
  stencil->SetStencil(dataToStencil->GetOutput());
  stencil->SetBackgroundValue(value);

  if (reverse == true){
    stencil->ReverseStencilOn();
  } else {
    stencil->ReverseStencilOff();
  }


  vtkImageData *vtkimageOut = vtkImageData::New();
  vtkimageOut = stencil->GetOutput();
  vtkimageOut->Modified();
  vtkimageOut->Update();

  // Retrieve the output in irtk canonical format.
  n    = imageCanonical.GetNumberOfVoxels();
  ptr1 = imageCanonical.GetPointerToVoxels();
  ptr2 = (short *)vtkimageOut->GetScalarPointer();
  for (i = 0; i < n; i++){
    *ptr1 = *ptr2;
    ptr1++;
    ptr2++;
  }

  // Prepare output image by clearing it.
  n    = image.GetNumberOfVoxels();
  ptr1 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++){
    *ptr1 = 0;
    ptr1++;
  }

  // Paste the canonical data.
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
        	val = imageCanonical.Get(i-i1, j-j1, k-k1);
        	image.Put(i, j, k, val);
        }
      }
    }

  image.Write(output_name);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}

#endif

/// BEGIN GRAVEYARD:

//   dataToStencil->SetInputConnection(stripper->GetOutputPort());
//   dataToStencil->Update();

//   stencil->Update();

//   for (i = 0; i < 6; ++i){
//     cout << bounds[i] << " ";
//   }
//   cout <<endl;

//   cout << i1 << " " << j1 << " " << k1 << " " << i2 << " " << j2 << " " << k2 << endl;

//   cout << "After transforming polys:" << endl;
//   polys->GetBounds(bounds);
//   for (int i = 0; i < 6; ++i){
//     cout << bounds[i] << " ";
//   }
//   cout << endl;

//   vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//   writer->SetInput(polys);
//   writer->SetFileName("newPolys.vtk");
//   writer->Write();

//   image3.Write("image3.nii.gz");
//   exit(0);

//   image.Write("bla.nii.gz");


//   cout << "vtkimage bounds:" << endl;
//   vtkimage->GetExtent(i1, i2, j1, j2, k1, k2);
//   cout << i1 << " " << i2 << " " << j1 << " " << j2 << " " << k1 << " " << k2 << endl;
//   cout << endl;

//   vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
//   triFilter->SetInput(polys);
//   triFilter->Update();


//   vtkStripper *stripper = vtkStripper::New();
//   stripper->SetInput(triFilter->GetOutput());
//   stripper->Update();

//   cout << "vtkimageOut bounds:" << endl;
//   vtkimageOut->GetExtent(i1, i2, j1, j2, k1, k2);
//   cout << i1 << " " << i2 << " " << j1 << " " << j2 << " " << k1 << " " << k2 << endl;
//   cout << endl;

//   double scalarRange[2];
//   vtkimageOut->GetScalarRange(scalarRange);
//   cout << "vtkimageOut range : " << scalarRange[0] << " " << scalarRange[1] << endl;
//   cout << "vtkimageOut ScalarTypeAsString : " << vtkimageOut->GetScalarTypeAsString() << endl;

//   cout << "1" << endl;
//   vtkImageStencilData *stencilData = vtkImageStencilData::New();

//   double ox, oy, oz, dx, dy, dz;
//   vtkimage->GetOrigin(ox, oy, oz);
//   stencilData->SetOrigin(ox, oy, oz);
//   vtkimage->GetSpacing(ox, oy, oz);
//   stencilData->SetSpacing(ox, oy, oz);
//   stencilData = dataToStencil->GetOutput();
//   stencilData->UpdateData();
//   cout << "2" << endl;

//   int iter;
//   int xmin, xmax, r1, r2;
//   /// / Initialize the VTK image
//   xmin = 0;
//   xmax = xdim - 1;

//   for (j = 0; j < ydim; ++j){
//     for (k = 0; k < zdim; ++k){
//       // Initialise iterator.
//       iter = 0;
//       stencilData->GetNextExtent(r1, r2, xmin, xmax, j, k, iter);
//       while (iter != 0){
//         for (i = r1; i <= r2; ++i){
//           image.Put(i, j, k, 1);
//         }
//         stencilData->GetNextExtent(r1, r2, xmin, xmax, j, k, iter);
//       }
//     }
//   }

/// END GRAVEYARD







