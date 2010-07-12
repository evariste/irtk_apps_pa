#if (defined HAS_VTK)

#include <irtkImage.h>

#include <nr.h>

#include <vtkFloatArray.h>
#include <vtkTriangle.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>


#include <vtkCurvatures.h>

char *input_name = NULL;
char *mask_name = NULL;


void usage()
{
  cerr << "Usage: polydatacurvatureindices [in] <options>" << endl;
  cerr << "" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
	int ok, i, j, invert = False;
  int noOfPoints;
  double K, H;
  vtkTriangle *triangle;
  vtkIdType *cells;
  unsigned short noOfCells;
  vtkIdList *ptIds;
  double v1[3], v2[3], v3[3];
  int nonTriangleFaces = 0;
  int unmaskedCount;
  double pMin = 0, pMax = 100;
  double kLo = FLT_MAX * -1.0;
  double hLo = FLT_MAX * -1.0;
  double kHi = FLT_MAX;
  double hHi = FLT_MAX;

  if (argc < 2){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;


  // Parse remaining arguments
  while (argc > 1){
    ok = False;
    if ((!ok) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-mask") == 0)){
    	argc--;
    	argv++;
    	mask_name = argv[1];
    	argc--;
    	argv++;
    	ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-pLo") == 0)){
    	argc--;
    	argv++;
    	pMin = atof(argv[1]);
    	argc--;
    	argv++;
    	ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-pHi") == 0)){
    	argc--;
    	argv++;
    	pMax = atof(argv[1]);
    	argc--;
    	argv++;
    	ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-opt") == 0)) {
      argc--;
      argv++;
      //	Do stuff and maybe
//      argc--;
//      argv++;
      ok = True;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData *surface;

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();

  surface = surface_reader->GetOutput();
  surface->Update();
  surface->BuildCells();
  surface->BuildLinks();

  noOfPoints = surface->GetNumberOfPoints();

  vtkCurvatures *curve_K = vtkCurvatures::New();
  curve_K->SetInput(surface);
  curve_K->SetCurvatureTypeToGaussian();

  vtkCurvatures *curve_H = vtkCurvatures::New();
  curve_H->SetInput(surface);
  curve_H->SetCurvatureTypeToMean();


  if (invert  == True){
    cout << "Inverting mean curvature flag." << endl;
    int invertFlag = curve_H->GetInvertMeanCurvature();
    curve_H->SetInvertMeanCurvature(1 - invertFlag);
  }

  // Retrieve the computed curvature values.
  curve_K->Update();
  curve_H->Update();

  // Deal with mask.
  int ind;
  int *mask = new int[noOfPoints];
  for (i = 0; i < noOfPoints; ++i){
  	mask[i] = 1;
  }

  if (mask_name != NULL){
    vtkFloatArray *mask_scalars = (vtkFloatArray*) surface->GetPointData()->GetArray(mask_name, ind);
    if (ind == -1 || mask_scalars == NULL){
      cerr << "Masking scalars unavailable with name " << mask_name << endl;
      exit(1);
    }
    if (mask_scalars->GetNumberOfComponents() > 1){
    	cerr << "Masking scalars " << mask_scalars->GetName() << " has more than one component." << endl;
    	exit(1);
    }

    for (i = 0; i < noOfPoints; ++i){
      if (mask_scalars->GetTuple1(i) > 0){
      	continue;
      }
      mask[i] = 0;
    }
  }



  vtkFloatArray *scalars_K = vtkFloatArray::New();
  vtkFloatArray *scalars_H = vtkFloatArray::New();

  scalars_K =
  		(vtkFloatArray*) curve_K->GetOutput()->GetPointData()->GetScalars("Gauss_Curvature");
  scalars_H =
  		(vtkFloatArray*) curve_H->GetOutput()->GetPointData()->GetScalars("Mean_Curvature");


  // Calculate areas.
  double *area = new double[noOfPoints];
  unmaskedCount = 0;

  for (i = 0; i < noOfPoints; ++i){
  	area[i] = 0.0;

  	if (mask[i] <= 0){
  		continue;
  	}

  	++unmaskedCount;

    surface->GetPointCells(i, noOfCells, cells);

    if ( cells == NULL )
      continue;

    for (j = 0; j < noOfCells; ++j){
      triangle = vtkTriangle::SafeDownCast(surface->GetCell(cells[j]));

      if ( triangle != NULL ){
        ptIds = triangle->GetPointIds();

        surface->GetPoint(ptIds->GetId(0), v1);
        surface->GetPoint(ptIds->GetId(1), v2);
        surface->GetPoint(ptIds->GetId(2), v3);

        area[i] += vtkTriangle::TriangleArea(v1, v2, v3) / 3.0;

      } else {
      	++nonTriangleFaces;
      }
    }
  }

  // Outlier control.
  float *kCentileData = new float[1 + unmaskedCount];
  float *hCentileData = new float[1 + unmaskedCount];

  unmaskedCount = 0;
  for (i = 0; i < noOfPoints; ++i){

  	if (mask <= 0)
  		continue;

  	kCentileData[unmaskedCount + 1] = scalars_K->GetTuple1(i);
  	hCentileData[unmaskedCount + 1] = scalars_H->GetTuple1(i);
  	++unmaskedCount;
  }

  sort(unmaskedCount, kCentileData);
  sort(unmaskedCount, hCentileData);

  ind = 1 + (int) round((double) pMin * (unmaskedCount-1) / 100.0);
  kLo = kCentileData[ind];
  hLo = hCentileData[ind];

  ind = 1 + (int) round((double) pMax * (unmaskedCount-1) / 100.0);
  kHi = kCentileData[ind];
  hHi = hCentileData[ind];



  // Now calculate indices.
  double muH = 0.0;
  double muK = 0.0;
  double MLN = 0.0;
  double ICI = 0.0;
  double GLN = 0.0;
  double ECI = 0.0;
  double totalArea = 0.0;
  double K2, H2;

  for (i = 0; i < noOfPoints; ++i){
    K = scalars_K->GetTuple1(i);
    H = scalars_H->GetTuple1(i);

    if (mask <= 0 || K < kLo || K > kHi || H < hLo || H > hHi)
    	continue;

    H2 = H*H;
    K2 = K*K;

    muK += K * area[i];
    muH += H * area[i];

    MLN += H2 * area[i];

    if (K > 0){
    	ICI += K * area[i];
    }

    GLN += K2 * area[i];

    if (H2 > K)
    	ECI += 4 * H * sqrt(H2 - K) * area[i];

    totalArea += area[i];
  }

  // Finishing off some of the calculations.
  muH = muH / totalArea;
  muK = muK / totalArea;

  MLN = sqrt(MLN / 4.0 / M_PI);

  GLN = sqrt(totalArea * GLN);


  cout << muH << endl;
  cout << muK << endl;
  cout << MLN << endl;
  cout << ICI << endl;
  cout << GLN << endl;
  cout << ECI << endl;
  cout << totalArea << endl;
  cout << kLo << " " << kHi << endl;
  cout << hLo << " " << hHi << endl;


  delete [] area;
  delete [] kCentileData;
  delete [] hCentileData;

  return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
