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
char *output_label_image_name = NULL;
char *output_scalars_name = NULL;

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
  cerr << "Assign labels to the vertices of surfaceIn.  The value assigned" << endl;
  cerr << "to a vertex is the label of the nearest voxel in labelImage." << endl;
  cerr << "" << endl;
  cerr << "Write the result to surfaceOut." << endl;
  cerr << "Options: " << endl;
  cerr << "" << endl;
  cerr << "  name [name]                 : Name of output scalar array in surface." << endl;
  cerr << "  -writeDilatedLabels [name]  : Write out image of dilated labels." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, noOfPoints;
  bool ok;
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
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-writeDilatedLabels") == 0)){
      argc--;
      argv++;
      output_label_image_name  = argv[1];
      argc--;
      argv++;
      // do stuff and maybe argv++ etc.
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
      argc--;
      argv++;
      output_scalars_name  = argv[1];
      argc--;
      argv++;
      // do stuff and maybe argv++ etc.
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Reading image ..." << endl;
  irtkGreyImage *labelImage = new irtkGreyImage(input_image_name);

  vtkPolyData *surface = vtkPolyData::New();

  // Read surface
  cerr << "Reading surface ... " << endl;
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_surface_name);
  reader->Update();
  surface = reader->GetOutput();

  noOfPoints= surface->GetNumberOfPoints();

  vtkIntArray *labelsOut = vtkIntArray::New();
  labelsOut->SetNumberOfComponents(1);
  labelsOut->SetNumberOfTuples(noOfPoints);


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


  if (output_label_image_name != NULL){
    dilatedLabels->Write(output_label_image_name);
  }


  cerr << "Assigning scalars using dilated labels ... " << endl;
  irtkImageFunction *interp = NULL;
  interp = new irtkNearestNeighborInterpolateImageFunction;
  interp->SetInput(dilatedLabels);
  interp->Initialize();
  zeroCount = 0;
  for (i = 0; i < noOfPoints; ++i){
    surface->GetPoint(i, pt);
    dilatedLabels->WorldToImage(pt[0], pt[1], pt[2]);

    val = (int) round(interp->Evaluate(pt[0], pt[1], pt[2]));
    labelsOut->SetTuple1(i,val);

  }


  cerr << "Updating surface ... " << endl;
  labelsOut->Modified();

  if (output_scalars_name != NULL){
    labelsOut->SetName(output_scalars_name);
  } else {
    labelsOut->SetName("Labels");
  }

  surface->GetPointData()->AddArray(labelsOut);
  if (output_scalars_name != NULL){
    surface->GetPointData()->SetActiveScalars(output_scalars_name);
  } else {
    surface->GetPointData()->SetActiveScalars("Labels");
  }
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

