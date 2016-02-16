#ifdef HAS_VTK


#include <irtkImage.h>


#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencil.h>
#include <vtkImageData.h>

#include <irtkTransformation.h>


char *input_polys_name = NULL, *output_name = NULL;
char *template_image_name = NULL;

void usage()
{
  cerr << "polydata2maskimage [polydata] [template image] [output image] <options>" << endl;
  cerr << "-value val : value to put into the structure (default 255)." << endl;
  cerr << "-reverse   : toggle whether the structure or the background gets the value." << endl;
  cerr << "" << endl;
  cerr << "Based on code at www.cmake.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataToImageData" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int i;
  irtkGreyPixel *ptr1;
  double bounds[6];
  short outValue = 1000;
  int reverse = false;
  double xspacing, yspacing, zspacing;
  double spacing[3];
  int dim[3];
  double minSpacing;


  if (argc < 4){
    usage();
  }

  // Parse parameters
  input_polys_name = argv[1];
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
      outValue = atoi(argv[1]);
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

  ////////////////////////////////////////

  // Read template image
  irtkGreyImage templateImg;
  templateImg.Read(template_image_name);


  // Read surface
  vtkPolyDataReader *surfReader = vtkPolyDataReader::New();
  surfReader->SetFileName(input_polys_name);
  surfReader->Modified();
  surfReader->Update();

  vtkTriangleFilter *triFilter = vtkTriangleFilter::New();
  triFilter->SetInputData(surfReader->GetOutput());
  triFilter->Update();

  vtkPolyData *polys = triFilter->GetOutput();


  // Image to be stenciled will have a smaller resolution than
  // the template and will cover the bounds of the mesh.

  vtkSmartPointer<vtkImageData> whiteImage =
    vtkSmartPointer<vtkImageData>::New();

  polys->GetBounds(bounds);

  char s[] = "xyz";
  cout << "Polydata bounds: " << endl;
  for (int i = 0; i < 3; i++){
    cout << " " <<  s[i] << ": (" << bounds[i * 2] << ", " << bounds[i * 2+1] << ")" << endl;
  }

  // Use a spacing smaller than the minimum voxel spacing in the original image.
  templateImg.GetPixelSize(&xspacing, &yspacing, &zspacing);

  minSpacing = min(xspacing, yspacing);
  minSpacing = min(minSpacing, zspacing);

  spacing[0] = 0.3 * minSpacing;
  spacing[1] = 0.3 * minSpacing;
  spacing[2] = 0.3 * minSpacing;
  whiteImage->SetSpacing(spacing);

  // compute dimensions
  for (int i = 0; i < 3; i++){
    dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
  }

  whiteImage->SetDimensions(dim);
  whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

  // vtk origin is the voxel at image indices 0,0,0
  double origVTK[3];
  origVTK[0] = bounds[0] + spacing[0] / 2;
  origVTK[1] = bounds[2] + spacing[1] / 2;
  origVTK[2] = bounds[4] + spacing[2] / 2;

  whiteImage->SetOrigin(origVTK);

  whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  // fill the image with foreground voxels:
  unsigned char fgVal = 255;
  unsigned char bgVal = 0;

  vtkIdType count = whiteImage->GetNumberOfPoints();

  for (vtkIdType i = 0; i < count; ++i){
    whiteImage->GetPointData()->GetScalars()->SetTuple1(i, fgVal);
  }

  // polygonal data --> image stencil:
  vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();

  pol2stenc->SetInputData(polys);

  pol2stenc->SetOutputOrigin(origVTK);
  pol2stenc->SetOutputSpacing(spacing);
  pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  pol2stenc->Update();

  // Cut the corresponding white image and set the background:
  vtkSmartPointer<vtkImageStencil> imgstenc =
    vtkSmartPointer<vtkImageStencil>::New();

  imgstenc->SetInputData(whiteImage);
  imgstenc->SetStencilData(pol2stenc->GetOutput());

  if (reverse == true){
    imgstenc->ReverseStencilOn();
  } else {
    imgstenc->ReverseStencilOff();
  }

  imgstenc->SetBackgroundValue(bgVal);
  imgstenc->Update();

  // Image to collect result at the higher resolution
  vtkSmartPointer<vtkImageData> imgHiResVTK = vtkImageData::New();

  imgHiResVTK = imgstenc->GetOutput();
  imgHiResVTK->Modified();


  // Attributes for IRTK nifti image.
  irtkImageAttributes attr;
  attr._t = 1;
  attr._x = dim[0];
  attr._y = dim[1];
  attr._z = dim[2];

  attr._dt = 1;
  attr._dx = spacing[0];
  attr._dy = spacing[1];
  attr._dz = spacing[2];

  attr._xorigin = origVTK[0] + (dim[0] - 1) * spacing[0] / 2.0;
  attr._yorigin = origVTK[1] + (dim[1] - 1) * spacing[1] / 2.0;
  attr._zorigin = origVTK[2] + (dim[2] - 1) * spacing[2] / 2.0;

  // Retrieve data from stencil output to an IRTK image.
  irtkGreyImage *imgHiRes = new irtkGreyImage(attr);

  vtkTypeUInt8 *ptr2;

  ptr1 = imgHiRes->GetPointerToVoxels();
  ptr2 = (vtkTypeUInt8 *)imgHiResVTK->GetScalarPointer();

  for (i = 0; i < imgHiRes->GetNumberOfVoxels(); i++){

    if (*ptr2 > 0)
      *ptr1 = outValue;

    ptr1++;
    ptr2++;
  }



  // Transform the high resolution version onto the required template.
  // Use linear interpolator.
  irtkImageTransformation *imagetransformation =
      new irtkImageTransformation;
  irtkTransformation *transformation = new irtkRigidTransformation;
  irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;

  imagetransformation->SetInput(imgHiRes, transformation);
  imagetransformation->SetOutput(&templateImg);
  imagetransformation->PutInterpolator(interpolator);
  imagetransformation->Run();

  // Threshold at 50%
  int halfVal = outValue / 2;

  ptr1 = templateImg.GetPointerToVoxels();
  for (i = 0; i < templateImg.GetNumberOfVoxels(); i++){
    if (*ptr1 >= halfVal){
      *ptr1 = outValue;
    } else {
      *ptr1 = 0;
    }

    ptr1++;
  }


  templateImg.Write(output_name);



}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}

#endif






