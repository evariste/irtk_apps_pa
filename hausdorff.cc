#ifdef HAS_VTK

///////////////////////////////////////////////////////////////
// Find a Hausdorff type distance betweeen a pair of surfaces.
// 

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <irtkLocator.h>
#include <nr.h>

char *surface_1_name = NULL;
char *surface_2_name = NULL;

void usage(){
  cerr << "\n   hausdorff [surface_1] [surface_2] <-options>\n" << endl; 
  cerr << "   Find a symmetric Hausdorff type distance betweeen a pair of surfaces." << endl;
  cerr << "   First, loop over all points in surface_1 and find the closest point in surface_2." << endl;
  cerr << "   The distances are ordered and the required percentile is found. " << endl;
  cerr << "   The process is repeated with the roles of surface_1 and surface_2 reversed. " << endl;
  cerr << "   The average of the values found is returned. " << endl;
  cerr << "   Options: " << endl;
  cerr << "     <-locator>   Locator : 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
  cerr << "     <-percentile>        : The percentile to use when finding the distance from one surface" << endl;
  cerr << "                            to the other,  default = 100, i.e. maximum." << endl;
  cerr << "   (The conventional Hausdorff distance is taken as max{H(1, 2), H(2, 1)}" << endl;
  cerr << "      Where H(A, B) =  max   {  min   dist(p, q) }" << endl;
  cerr << "                      p in A   q in B" << endl;
  cerr << "    i.e. to get the conventional Hausdorff distance leave percentile to default value" << endl;
  cerr << "    and take max of H(1, 2) and H(2, 1) )" << endl << endl;;
  exit(1);

}

int main(int argc, char **argv){

  int i, ok, locatorType, percentile = 100;

  if (argc < 3){
    usage();
  }

  surface_1_name = argv[1];
  argv++;
  argc--;
  surface_2_name = argv[1];
  argv++;
  argc--;

  locatorType = 1;

  while (argc > 1){
    ok = False;

    if ((ok == False) && (strcmp(argv[1], "-locator") == 0)){
      argc--;	
      argv++;
      locatorType = atoi(argv[1]);
      argc--;	
      argv++;	
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-percentile") == 0)){
      argc--;	
      argv++;
      percentile = atoi(argv[1]);
      argc--;	
      argv++;	
      ok = True;
    }

    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  if (percentile < 0 || percentile > 100){
    cerr << "Invalid value for percentile given (" << percentile << ")." << endl;
    exit(1);
  }

  cout << "Reading data ... "; cout.flush();
  vtkPolyDataReader *reader1 = vtkPolyDataReader::New();
  reader1->SetFileName(surface_1_name);
  reader1->Update();

  vtkPolyData* surface_1 = vtkPolyData::New();
  surface_1 = reader1->GetOutput();

  vtkPolyDataReader *reader2 = vtkPolyDataReader::New();
  reader2->SetFileName(surface_2_name);
  reader2->Update();

  vtkPolyData* surface_2 = vtkPolyData::New();
  surface_2 = reader2->GetOutput();
  cout << "done" << endl;

  // Create locators
  irtkLocator *locator_1 = new irtkLocator;
  locator_1->SelectLocatorType(locatorType);
  locator_1->SetDataSet(surface_1);

  irtkLocator *locator_2 = new irtkLocator;
  locator_2->SelectLocatorType(locatorType);
  locator_2->SetDataSet(surface_2);

  int n_1 = surface_1->GetNumberOfPoints();
  int n_2 = surface_2->GetNumberOfPoints();

  double pt_a[3], pt_b[3];

  float *dist_1 = new float[1 + n_1];
  float *dist_2 = new float[1 + n_2];

  for (i = 0; i < n_1; i++){
    surface_1->GetPoints()->GetPoint (i, pt_a);

    pt_b[0] = pt_a[0];
    pt_b[1] = pt_a[1];
    pt_b[2] = pt_a[2];

    (void) locator_2->FindClosestPoint (pt_b);

    dist_1[1 + i] = sqrt((pt_a[0] - pt_b[0]) * (pt_a[0] - pt_b[0]) +
                         (pt_a[1] - pt_b[1]) * (pt_a[1] - pt_b[1]) +
                         (pt_a[2] - pt_b[2]) * (pt_a[2] - pt_b[2]));

  }

  for (i = 0; i < n_2; i++){
    surface_2->GetPoints()->GetPoint (i, pt_a);

    pt_b[0] = pt_a[0];
    pt_b[1] = pt_a[1];
    pt_b[2] = pt_a[2];

    (void) locator_1->FindClosestPoint (pt_b);

    dist_2[1 + i] = sqrt((pt_a[0] - pt_b[0]) * (pt_a[0] - pt_b[0]) +
                         (pt_a[1] - pt_b[1]) * (pt_a[1] - pt_b[1]) +
                         (pt_a[2] - pt_b[2]) * (pt_a[2] - pt_b[2]));

  }

  sort(n_1, dist_1);
  sort(n_2, dist_2);

  int index_1 = 1 + (int) round( (double) percentile * (n_1 - 1) / 100.0);
  int index_2 = 1 + (int) round( (double) percentile * (n_2 - 1)/ 100.0);

  double val_1 = dist_1[index_1];
  double val_2 = dist_2[index_2];

  cout << "  percentile " << percentile << endl;
  cout << "     H(1, 2) " << val_1 << endl;
  cout << "     H(2, 1) " << val_2 << endl;
  cout << "     Average " << 0.5 * (val_1 + val_2) << endl;

  delete [] dist_1;
  delete [] dist_2;

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
