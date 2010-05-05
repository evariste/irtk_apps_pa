#if (defined HAS_VTK)

#ifdef WIN32
#include <map>
#else
#include <map>
#endif

#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkEuclideanDistanceTransform.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

char *input_image_name = NULL, *output_name = NULL;
char *input_surface_name = NULL;

typedef map<short, long> countMap;

map<short, long>::iterator iter;

// Returns a floating point 0/1 image to show a label.  The image needs to
// be floating point so that it can later be used in a distance map filter.
irtkRealImage *getBinaryLabelImage(irtkGreyImage *labelImage, short label){

  int i, noOfVoxels;
  irtkRealPixel *ptr2mask;
  irtkGreyPixel *ptr2labels;
  irtkRealImage *mask = new irtkRealImage(*labelImage);
  ptr2mask = mask->GetPointerToVoxels();
  ptr2labels = labelImage->GetPointerToVoxels();

  noOfVoxels = labelImage->GetNumberOfVoxels();

  for (i = 0; i < noOfVoxels; ++i, ++ptr2mask, ++ptr2labels){
    if (*ptr2labels == label){
      *ptr2mask = 1;
    } else {
      *ptr2mask = 0;
    }
  }
  return mask;
}

void usage()
{
  cerr << "Usage:  polydataassignlabels [labelImage] [surfaceIn] [surfaceOut] <options>" << endl;
  cerr << "" << endl;
  cerr << "Assign scalars to the vertices of surfaceIn.  The scalar assigned" << endl;
  cerr << "is the label of the nearest voxel in labelImage to the vertex." << endl;
  cerr << "" << endl;
  cerr << "Write the result to surfaceOut." << endl;
  cerr << "Options: " << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, ok, noOfPoints;
  double surfaceBounds[6];
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double pt[3];
  int val, zeroCount;
  int noOfVoxels, noOfLabels, currLabel;

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
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-XXX") == 0)){
      argc--;
      argv++;
      // do stuff and maybe argv++ etc.
      ok = True;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Reading image ..." << endl;
  irtkGreyImage *labelImage = new irtkGreyImage(input_image_name);

  vtkPolyData *surface = vtkPolyData::New();

  cerr << "Reading surface ... " << endl;
  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_surface_name);
  reader->Update();
  surface = reader->GetOutput();

  noOfPoints= surface->GetNumberOfPoints();

  //cerr << "Making space for scalars ..., no. of points =  " << noOfPoints << endl;
  vtkIntArray *scalars = vtkIntArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(noOfPoints);

  //cerr << "Creating interpolator " << endl;
  irtkImageFunction *interp = NULL;
  interp = new irtkNearestNeighborInterpolateImageFunction;
  interp->SetInput(labelImage);
  interp->Initialize();

  //Check that surface does not go outside fov of label image.
  surface->ComputeBounds();
  surface->GetBounds(surfaceBounds);
  xmin = surfaceBounds[0];
  xmax = surfaceBounds[1];
  ymin = surfaceBounds[2];
  ymax = surfaceBounds[3];
  zmin = surfaceBounds[4];
  zmax = surfaceBounds[5];
  cerr << "Bounds of surface : ";
  cerr << "(" << xmin << ", " << ymin << ", " << zmin << ") and ";
  cerr << "(" << xmax << ", " << ymax << ", " << zmax << ")" << endl;

  labelImage->WorldToImage(xmin, ymin, zmin);
  labelImage->WorldToImage(xmax, ymax, zmax);
  cerr << "In image coords : ";
  cerr << "(" << xmin << ", " << ymin << ", " << zmin << ") and ";
  cerr << "(" << xmax << ", " << ymax << ", " << zmax << ")" << endl;
  if (xmin < -0.5 || xmax > labelImage->GetX()-0.5 ||
      ymin < -0.5 || ymax > labelImage->GetY()-0.5 ||
      zmin < -0.5 || zmax > labelImage->GetZ()-0.5){
        cerr << "Surface outside bounds of image." << endl;
        exit(1);
  }

  irtkRealImage *currLabelMask;
  irtkGreyImage *dilatedLabels;
  irtkRealImage *minDmap;
  irtkRealImage *currDmap;

  // Identify the number of distinct labels in the image.
  irtkRealPixel *ptr2dmap;
  irtkRealPixel *ptr2minDmap;
  irtkGreyPixel *ptr2label;
  irtkGreyPixel *ptr2dilatedLabel;

  countMap labelCount;

  // Count up different labels so we can identify the number of distinct labels.
  noOfVoxels = labelImage->GetNumberOfVoxels();
  ptr2label  = labelImage->GetPointerToVoxels();
  for (i = 0; i < noOfVoxels; ++i, ++ptr2label){
    if (*ptr2label > 0){
      ++labelCount[*ptr2label];
    }
  }

  cerr << "No. of Voxels          : " << noOfVoxels << endl;
  cerr << "Label Counts " << endl;
  for (iter = labelCount.begin(); iter != labelCount.end(); ++iter){
    cerr << iter->first << "\t" << iter->second << endl;
  }
  noOfLabels = labelCount.size();
  cerr << "No. of distinct labels : " << noOfLabels << endl;

  // Using the distance maps.
  minDmap  = new irtkRealImage(*labelImage);
  currDmap = new irtkRealImage(*labelImage);

  // Note that the dilated labels are initialised to the given label image.
  // I.e. the original labels are left alone and we seek to assign labels to
  // the zero voxels based on closest labeled voxels.
  dilatedLabels = new irtkGreyImage(*labelImage);

  // Initialise the minimum distance map.
  ptr2minDmap = minDmap->GetPointerToVoxels();
  for (i = 0; i < noOfVoxels; ++i, ++ptr2minDmap){
    *ptr2minDmap = FLT_MAX;
  }

  // Single distance transform filter for all labels.
  irtkEuclideanDistanceTransform<irtkRealPixel> *edt = NULL;
  edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
    (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);

  cerr << "Finding distance maps ..." << endl;
  cerr << "Current label : " << endl;
  for (iter = labelCount.begin(); iter != labelCount.end(); ++iter){
    currLabel = iter->first;
    cerr << "  " << currLabel << endl;

    // There is a new operator used for currLabelMask in the following function call.
    currLabelMask = getBinaryLabelImage(labelImage, currLabel);

    edt->SetInput(currLabelMask);
    edt->SetOutput(currDmap);
    edt->Run();

    ptr2minDmap      = minDmap->GetPointerToVoxels();
    ptr2dmap         = currDmap->GetPointerToVoxels();
    ptr2label        = labelImage->GetPointerToVoxels();
    ptr2dilatedLabel = dilatedLabels->GetPointerToVoxels();

    for (i = 0; i < noOfVoxels; ++i, ++ptr2minDmap, ++ptr2dmap, ++ptr2label, ++ptr2dilatedLabel){
      if (*ptr2label == 0 && *ptr2dmap < *ptr2minDmap){
        *ptr2minDmap = *ptr2dmap;
        *ptr2dilatedLabel = currLabel;
      }
    }

    // Tidy up.
    delete currLabelMask;
  }

  cerr << "Assigning scalars using dilated labels ... " << endl;
  interp->SetInput(dilatedLabels);
  interp->Initialize();
  zeroCount = 0;
  for (i = 0; i < noOfPoints; ++i){
    surface->GetPoint(i, pt);
    dilatedLabels->WorldToImage(pt[0], pt[1], pt[2]);

    val = (int) round(interp->Evaluate(pt[0], pt[1], pt[2]));
    scalars->SetTuple1(i,val);

    //if (val == 0){
    //  ++zeroCount;
    //}
  }

  //cerr << "Zero count : " << zeroCount << " = " << 100.0 * zeroCount / ((double) noOfPoints) << "%" << endl;
  cerr << "Updating surface ... " << endl;
  scalars->Modified();
  scalars->SetName("Labels");
  surface->GetPointData()->SetScalars(scalars);
  surface->Update();

  cerr << "Writing surface ... " << endl;
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

