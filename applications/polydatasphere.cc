#if (defined HAS_VTK)

#include <irtkImage.h>

#include <map>

#include <nr.h>
#include <time.h>

#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
//#include <vtkSphereSource.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

#ifndef M_PI_2
# define M_PI_2		1.57079632679489661923	// pi/2
#endif

typedef struct {
    int id0, id1;
} edge;

struct edgeCmp {
   bool operator()( const edge ea, const edge eb ) const {
     return ( ea.id0 < eb.id0 ) || (( ea.id0 == eb.id0 ) && ( ea.id1 < eb.id1));
   }
 };


typedef map<edge, int, edgeCmp> midpointMap;

map<edge, short>::iterator iter;



#define MAXVALS 100

char *input_name = NULL;
char *output_name = NULL;
char *inputScalarName = NULL;

double re_ylm(int l, int m, double theta, double phi);

double mag_ylm(int l, int m, double theta, double phi);
double plgndr2(int l, int m, double x);

int getEdgeCount(vtkIdList** adj, int noOfVertices);
void getAdjacency(vtkCellArray *faces, vtkIdList** adj);
void addAdjacency(vtkIdList**adj, vtkIdType a, vtkIdType b);
void deleteIdListArray(vtkIdList** &list, int size);
void copyIdListArray(vtkIdList** &src, vtkIdList** &dest, int noOfVerts);
void createIdListArray(vtkIdList** &list, int size);
void copyFaces(vtkCellArray *src, vtkCellArray *dest);
void copyPoints(vtkPoints *src, vtkPoints*dest);
void midpoint(double *x, double *y, double *m);

void get_ylm_values(vtkFloatArray* theta, vtkFloatArray* phi, vtkFloatArray* ylm, int l, int m, double &min, double &max);

vtkPolyData *getSphereApproximation();

void usage()
{
  cerr << " polydatabumps [output] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Generate a random set of bumps on a unit sphere using a set of" << endl;
  cerr << " spherical harmonics up to a chosen degree." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -l [val]  The maximum degree of spherical harmonics to be used." << endl;
  cerr << "           Higher l gives more complex surfaces but takes longer"   << endl;
  cerr << "           to compute." << endl;
  cerr << " -uniform  Use this flag to give equal weights to the basis functions " << endl;
  cerr << "           The default is to combine them with random weights. " << endl;
  cerr << " -r        The size of the bumps." << endl;

  exit(1);
}



int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  int i;
  bool ok, useRand;
  int noOfPoints;
  double x, y, z, r, pt[3];
  double a, b, bump_r;
  double val;
  double min, max, maxabs;
//  float t, p;
  time_t seconds;
  long ran2Seed;
  long ran2initialSeed;

  float *weight;
  double weightSum;

  seconds = time(NULL);
  ran2Seed = seconds;
  ran2initialSeed = -1 * ran2Seed;
  (void) ran2(&ran2initialSeed);

  useRand = true;

  a = 1.0f;
  b = 1.0f;
  bump_r = 0.5f;
  r = 1.0f;

  output_name  = argv[1];
  argc--;
  argv++;

  int l, l_max, m;
  l_max = 1;
  m = 1;


  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-l") == 0)){
      argc--;
      argv++;
      l_max = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-uniform") == 0)){
      argc--;
      argv++;
      useRand = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-r") == 0)){
      argc--;
      argv++;
      bump_r = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  ////////////////////

  vtkPolyData *output = vtkPolyData::New();

  output = getSphereApproximation();
  output->Update();

  ////////////////////
  noOfPoints = output->GetNumberOfPoints();

  cout << noOfPoints << endl;

  vtkFloatArray *phi = vtkFloatArray::New();
  phi->SetNumberOfComponents(1);
  phi->SetNumberOfTuples(noOfPoints);

  vtkFloatArray *theta = vtkFloatArray::New();
  theta->SetNumberOfComponents(1);
  theta->SetNumberOfTuples(noOfPoints);

  // Find azimuth (phi) and altitude (theta)
  for (i = 0; i < noOfPoints; ++i){
	  output->GetPoint(i, pt);
	  x = pt[0];
	  y = pt[1];
	  z = pt[2];

	  // Add pi to ensure range is in 0->2 pi
	  val = atan2f(y, x);

	  if (val < 0)
	    val += 2 * M_PI;

	  phi->SetTuple1(i, val);

	  r = sqrt(x*x + y*y);
	  // Adding pi/2 ensures range is from 0 to pi, angle measured from north
	  // pole.
	  val = atan2f(r, z);

	  if (val < 0)
	    val += M_PI_2;

    theta->SetTuple1(i, val);
  }


  phi->SetName("phi");
  output->GetPointData()->AddArray(phi);
  cout << "Added phi (azimuth) array." << endl;

  theta->SetName("theta");
  output->GetPointData()->AddArray(theta);
  cout << "Added theta (altitude) array." << endl;

  ///////////////////////////////////////////////////

  vtkFloatArray *ylm = vtkFloatArray::New();
  ylm->SetNumberOfComponents(1);
  ylm->SetNumberOfTuples(noOfPoints);

  vtkFloatArray *randFunc = vtkFloatArray::New();
  randFunc->SetNumberOfComponents(1);
  randFunc->SetNumberOfTuples(noOfPoints);
  randFunc->FillComponent(0, 0.0);

  int wIndex, wCount;

  wIndex = 0;
  for (l = 0; l <= l_max; ++l){
    for (m = 0; m <= l; ++m){
      ++wIndex;
    }
  }

  wCount = wIndex;
  cout << "l_max  : " << l_max << endl;
  cout << "wCount : " << wCount << endl;

//  weight = new float[l_max + 1];
  weight = new float[wCount];

  weightSum = 0.0;
  wIndex = 0;

  for (l = 0; l <= l_max; ++l){
    for (m = 0; m <= l; ++m){

      if (useRand == true){
        weight[wIndex] = ran2(&ran2Seed);
      } else {
        weight[wIndex] = 1.0;
      }

      weightSum += weight[wIndex];
      cerr << wIndex << " " << weight[wIndex] << " " << weightSum << endl;
      ++wIndex;
    }
  }

  wCount = wIndex;

//  sort(wCount, weight - 1);

  cerr << endl;

  wIndex = 0;
  for (l = 0; l <= l_max; ++l){
    for (m = 0; m <= l; ++m){

      get_ylm_values(theta, phi, ylm, l, m, min, max);

      // Add a weighted ylm to the function being accumulated.
      for (i = 0; i < noOfPoints; ++i){
        val = randFunc->GetTuple1(i);
        val += weight[wIndex] * ylm->GetTuple1(i);
        randFunc->SetTuple1(i, val);
      }

      ++wIndex;
    }
  }
  cerr << endl;

  min = FLT_MAX;
  max = -1.0 * FLT_MAX;
  for (i = 0; i < noOfPoints; ++i){
    val = randFunc->GetTuple1(i);
    val /= weightSum;

    if (val < min)
      min = val;
    if (val > max)
      max = val;

    randFunc->SetTuple1(i, val);
  }
  cout << "randFunc : " << min << "   " << max << endl;

  if (fabs(min) > max){
    maxabs = fabs(min);
  } else {
    maxabs = fabs(max);
  }

  cout << "max abs : " << maxabs << endl;

  /////////////////////////////////////////////////
  // Adding bumps.

  for (i = 0; i < noOfPoints; ++i){

    output->GetPoint(i, pt);
    x = pt[0];
    y = pt[1];
    z = pt[2];

    val = randFunc->GetTuple1(i);
//    val = -0.5 + (val - min) / (max - min);
//    val = val * val;
    val = fabs(val) / maxabs;
    val = sqrt(val);

    val = 1 + bump_r * val;
    // Update the function.
    randFunc->SetTuple1(i, val);

    x *= val;
    y *= val;
    z *= val;

    output->GetPoints()->SetPoint(i, x, y, z);

  }

  /////////////////////////////////////////////////

  randFunc->SetName("randFunc");
  output->GetPointData()->AddArray(randFunc);
  output->GetPointData()->SetActiveScalars("randFunc");

  cout << "Added randFunc array." << endl;

  /////////////////////////////////////////////////

  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();

  normalsFilter->SplittingOff();
  normalsFilter->SetInput(output);
  normalsFilter->Modified();
  normalsFilter->Update();

  /////////////////////////////////////////////////

  cerr << "Writing surface ... " << endl;
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(normalsFilter->GetOutput());
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  cerr << ".. done. " << endl;

  delete [] weight;

}


void get_ylm_values(vtkFloatArray* theta, vtkFloatArray* phi, vtkFloatArray* ylm, int l, int m, double &min, double &max)
{
  int i, noOfPoints;
  double t, p, val;

  noOfPoints = theta->GetNumberOfTuples();

  min = FLT_MAX;
  max = -1.0 * FLT_MAX;

  for (i = 0; i < noOfPoints; ++i){
    p = phi->GetTuple1(i);
    t = theta->GetTuple1(i);

    val = re_ylm(l, m, t, p);

    ylm->SetTuple1(i, val);

//    if (val < min)
//      min = val;
//    if (val > max)
//      max = val;
//

//  cout << "ylm : " << min << "   " << max << endl;
//  if (fabs(min) > max){
//    max = fabs(min);
  }

//  cout << "Overall max : " << max << endl;
//  cout << "--------------------" << endl;

}


double re_ylm(int l, int m, double theta, double phi)
{

  // Find degree l, order m, Legendre polynomial value at cos theta: plm(cos theta)
  int i;
  double val;

  double costheta;
  double plm_costheta;
  double cos_m_phi;

  val = 1;

  if (m < 0){
    // y_l_-m = -1^m conj(y_l_m)
    // we only need the real part so can ignore the conjugation part and just
    // start off our accumulated value with -1^m

    // So flip the sign of m
    m = -1 * m;

    // -1^m
    if (m % 2 == 0){
      val = 1;
    } else {
      val = -1;
    }
  }

  if (m > l){
    return 0;
  }

  // Now we are sure that 0 le m le l

  costheta = cos(theta);
//  plm_costheta = plgndr(l, m, costheta);
  plm_costheta = plgndr2(l, m, costheta);

//  if (plm_costheta > FLT_MAX){
//    cerr << "BIG " << "" << endl;
//    cerr << l << " " << m << " " << theta << " " << plm_costheta << endl;
//  }


  // Multiply by e^(i m phi) = cos(m phi) + i sin(m phi)
  // but ylm is complex valued, we only return the real part.

  cos_m_phi = cos(double (m *  phi));

  val *= plm_costheta;
  val *= cos_m_phi;

  // Multiply by the appropriate constants.
  val *= sqrt(2.0 * l + 1);
  val /= sqrt(4 * M_PI);

  int start, end;
  start = l - m + 1;
  end = l + m;

  for (i = start; i < end; ++i){
    val /= sqrt((double) i);
  }

  return val;
}


double mag_ylm(int l, int m, double theta, double phi){

  // Find degree l, order m, Legendre polynomial value at cos theta: plm(cos theta)
  int i;
  double val;

  double plm_costheta;

  val = 1;

  if (m < 0){
    // y_l_-m = -1^m conj(y_l_m)
    // we only need the magnitude so can ignore the conjugation part and can even
    // ignore the -1^m part

    // So just flip the sign of m
    m = -1 * m;
  }

  if (m > l){
    return 0;
  }

  // Now we are sure that 0 le m le l

//  plm_costheta = plgndr(l, m, cos(theta));
  plm_costheta = plgndr2(l, m, cos(theta));
  val *= plm_costheta;

  // Multiply by |e^(i m phi)| = 1 , i.e. no need to do anything.

  // Multiply by the appropriate constants.
  val *= sqrt(2.0 * l + 1);
  val /= sqrt(4 * M_PI);

  int start, end;
  start = l - m + 1;
  end = l + m;

  for (i = start; i < end; ++i){
    val /= sqrt((double) i);
  }

  return val;
}

double plgndr2(int l, int m, double x)
{
  void nrerror(char error_text[]);
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;

  pll = 0.0;

  if (m < 0 || m > l || fabs(x) > 1.0){
    cerr << "Bad arguments in routine plgndr2" << endl;
    exit(1);
  }
  pmm=1.0;
  if (m > 0) {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=m+2;ll<=l;ll++) {
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;
      }
      return pll;
    }
  }
}


vtkPolyData *getSphereApproximation(){

  int i, j, k;
  double a[3];
  double v0[3], v1[3], v2[3];
  vtkIdType npts = 0;
  vtkIdType id0, id1, id2, ida, idb, idc;
  vtkIdType *pts = NULL;
  int vertCounter;
  int noOfVertsOld = 0, noOfEdgesOld, noOfFacesOld;
  int noOfVertsNew, noOfEdgesNew, noOfFacesNew;
  midpointMap midptMap;
  edge e;
  int solidIndex = VTK_SOLID_ICOSAHEDRON;
  int levels;
  int maxPts, estSize;

  levels = 6;

  vtkPlatonicSolidSource *s = vtkPlatonicSolidSource::New();

  vtkPolyData *input = vtkPolyData::New();

  solidIndex = VTK_SOLID_ICOSAHEDRON;

  s->SetSolidType(solidIndex);

  input = s->GetOutput();
  input->Update();

  vtkPoints *oldPoints = vtkPoints::New();
  vtkCellArray* oldFaces = vtkCellArray::New();
  vtkIdList** oldAdj = NULL;

  vtkPoints *newPoints = vtkPoints::New();
  vtkCellArray* newFaces = vtkCellArray::New();
  vtkIdList** newAdj = NULL;

  copyPoints(input->GetPoints(), newPoints);

  copyFaces(input->GetPolys(), newFaces);


  noOfVertsNew = newPoints->GetNumberOfPoints();
  noOfFacesNew = newFaces->GetNumberOfCells();

  createIdListArray(newAdj, noOfVertsNew);

  // Convert the adjacency info. implicit in the faces into an adjacency list.
  getAdjacency(newFaces, newAdj);


  for (k = 0; k < levels; k++){

    // Delete old adjacency if necessary
    if (oldAdj != NULL){
      deleteIdListArray(oldAdj, noOfVertsOld);
    }

    // Copy the new points into the old (current) points.
    copyPoints(newPoints, oldPoints);
    noOfVertsOld = oldPoints->GetNumberOfPoints();

    // Same for faces.
    copyFaces(newFaces, oldFaces);
    noOfFacesOld = oldFaces->GetNumberOfCells();

    // Same for adjacency (edge) info.
    createIdListArray(oldAdj, noOfVertsOld);
    copyIdListArray(newAdj, oldAdj, noOfVertsOld);
    noOfEdgesOld = getEdgeCount(oldAdj, noOfVertsOld);
    deleteIdListArray(newAdj, noOfVertsOld);

    // How many vertices, edges and faces after the current iteration?
    // V' = V + E
    // F' = 4F
    // E' = 2E + 3F

//    cerr << "old:  v " << noOfVertsOld << " e " << noOfEdgesOld << " f " << noOfFacesOld << endl;

    noOfVertsNew = noOfVertsOld + noOfEdgesOld;
    noOfFacesNew = 4 * noOfFacesOld;
    noOfEdgesNew = 2 * noOfEdgesOld + 3 * noOfFacesOld;

    // Allocate space for new vertices edges faces
    newPoints->SetNumberOfPoints(noOfVertsNew);

    maxPts = oldFaces->GetMaxCellSize();
    newFaces->Initialize();
    estSize = newFaces->EstimateSize(noOfFacesNew, maxPts);
    newFaces->Allocate(estSize, 0);
    newFaces->Reset();

    createIdListArray(newAdj, noOfVertsNew);

    // Copy over original vertices, these will retain the same indices in the vertex array.
    for (i = 0; i < noOfVertsOld; ++i){
      oldPoints->GetPoint(i, a);
      newPoints->InsertPoint(i, a);
    }
    // Generate all the new vertices, one for each of the old edges.
    vertCounter = noOfVertsOld;
    for (i = 0; i < noOfVertsOld; ++i){
      id0 = i;
      newPoints->GetPoint(id0, v0);
      npts = oldAdj[id0]->GetNumberOfIds();
      // Loop over list of neighbours.
      for (j = 0; j < npts; ++j){
        id1 = oldAdj[i]->GetId(j);
        newPoints->GetPoint(id1, v1);
        midpoint(v0, v1, a);
        vtkMath::Normalize(a);
        newPoints->InsertPoint(vertCounter, a);

        // We know that id0 should be < id1
        e.id0 = id0;
        e.id1 = id1;
        midptMap[e] = vertCounter;
        ++vertCounter;
      }
    }

    // Loop over old faces, subdividing each.
    oldFaces->InitTraversal();
    newFaces->InitTraversal();

    for(i = 0; i < noOfFacesOld; ++i){
      oldFaces->GetNextCell(npts, pts);

      if (npts != 3){
        cerr << "Expecting triangles only!" << endl;
        exit(1);
      }

      id0 = pts[0];
      id1 = pts[1];
      id2 = pts[2];

      newPoints->GetPoint(id0, v0);
      newPoints->GetPoint(id1, v1);
      newPoints->GetPoint(id2, v2);

      /*
      Subdivide each triangle in the old approximation and normalize
      the new points thus generated to lie on the surface of the unit
      sphere.
      Each input triangle with vertices labelled [0,1,2] as shown
      below will be turned into four new triangles:

      Make new points
                      a = (0+2)/2
                      b = (0+1)/2
                      c = (1+2)/2
                   1
                   /\   Normalize a, b, c
                  /  \
                b/____\ c   Construct new triangles
                /\    /\        [0,b,a]
               /  \  /  \       [b,1,c]
              /____\/____\      [a,b,c]
             0    a       2     [a,c,2]
      */

      // Look up the ids of the midpoints.
      e.id0 = min(id0, id2);
      e.id1 = max(id0, id2);
      ida = midptMap[e];

      e.id0 = min(id0, id1);
      e.id1 = max(id0, id1);
      idb = midptMap[e];

      e.id0 = min(id1, id2);
      e.id1 = max(id1, id2);
      idc = midptMap[e];

      // Add the new adjacencies.
      addAdjacency(newAdj, id0, ida);
      addAdjacency(newAdj, id2, ida);

      addAdjacency(newAdj, id0, idb);
      addAdjacency(newAdj, id1, idb);

      addAdjacency(newAdj, id1, idc);
      addAdjacency(newAdj, id2, idc);

      addAdjacency(newAdj, ida, idb);
      addAdjacency(newAdj, ida, idc);
      addAdjacency(newAdj, idb, idc);

      // Add the new faces.
      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(id0);
      newFaces->InsertCellPoint(idb);
      newFaces->InsertCellPoint(ida);
      newFaces->UpdateCellCount(3);

      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(idb);
      newFaces->InsertCellPoint(id1);
      newFaces->InsertCellPoint(idc);
      newFaces->UpdateCellCount(3);

      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(ida);
      newFaces->InsertCellPoint(idb);
      newFaces->InsertCellPoint(idc);
      newFaces->UpdateCellCount(3);

      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(ida);
      newFaces->InsertCellPoint(idc);
      newFaces->InsertCellPoint(id2);
      newFaces->UpdateCellCount(3);

    }

    newFaces->Squeeze();
//    cerr << "new faces : " << newFaces->GetNumberOfCells() << endl;

  }

  vtkPolyData *output = vtkPolyData::New();
  output->SetPoints(newPoints);
  output->SetPolys(newFaces);
  output->Update();

  return output;

}

void midpoint(double *x, double *y, double *m){
  m[0] = 0.5 * (x[0] + y[0]);
  m[1] = 0.5 * (x[1] + y[1]);
  m[2] = 0.5 * (x[2] + y[2]);
}

void copyPoints(vtkPoints *src, vtkPoints*dest){
  int i, noOfPts;
  double p[3];

  noOfPts = src->GetNumberOfPoints();
  dest->SetNumberOfPoints(noOfPts);

  for (i = 0; i < noOfPts; ++i){
    src->GetPoint(i, p);
    dest->InsertPoint(i, p);
  }
//  cerr << "Copied " << noOfPts << " points" << endl;
}

void copyFaces(vtkCellArray *src, vtkCellArray *dest){
  int i, j, noOfCells;
  vtkIdType npts = 0;
  vtkIdType *ptIds = NULL;
  noOfCells = src->GetNumberOfCells();
  dest->SetNumberOfCells(noOfCells);

  int estSize;
  int maxPts;
  maxPts = src->GetMaxCellSize();

  dest->Initialize();
  estSize = dest->EstimateSize(noOfCells, maxPts);
  dest->Allocate(estSize, 0);
  dest->Reset();

  dest->InitTraversal();
  src->InitTraversal();
  for (i = 0; i < noOfCells; ++i){
    src->GetNextCell(npts, ptIds);
    dest->InsertNextCell(npts);
    for (j = 0; j < npts; ++j){
      dest->InsertCellPoint(ptIds[j]);
    }
    dest->UpdateCellCount(npts);
  }
  dest->Squeeze();
}

void createIdListArray(vtkIdList** &list, int size){
  int i;
  list = new vtkIdList*[size];
  for (i = 0; i < size; ++i){
    list[i] = vtkIdList::New();
  }
}

void copyIdListArray(vtkIdList** &src, vtkIdList** &dest, int noOfVerts){
  int i, j, nPts;
  // Assume both lists initialised to same size.
  for (i = 0; i < noOfVerts; ++i){
    nPts = src[i]->GetNumberOfIds();
    for (j = 0; j < nPts; ++j){
      dest[i]->InsertUniqueId(src[i]->GetId(j));
    }
  }
}

void deleteIdListArray(vtkIdList** &list, int size){
  int i;
  for (i = 0; i < size; ++i){
    list[i]->Delete();
  }
  delete [] list;
  list = NULL;
}

void addAdjacency(vtkIdList**adj, vtkIdType a, vtkIdType b){
  int temp;
  if (a == b)
    return;
  if (a > b){
    temp = a;
    a = b;
    b = temp;
  }
  adj[a]->InsertUniqueId(b);
}

void getAdjacency(vtkCellArray *faces, vtkIdList** adj){

  int i, j, noOfFaces;
  vtkIdType ida, idb;

  vtkIdType npts = 0;
  vtkIdType *pts = NULL;

  noOfFaces = faces->GetNumberOfCells();

  faces->InitTraversal();
  for (i = 0; i < noOfFaces; ++i){

    faces->GetNextCell(npts, pts);

    ida = pts[0];
    idb = pts[npts - 1];

    if (ida < idb){
      adj[ida]->InsertUniqueId(idb);
    } else {
      adj[idb]->InsertUniqueId(ida);
    }
    for (j = 0; j < npts - 1; ++j){
      ida = pts[j];
      idb = pts[j + 1];
      if (ida < idb){
        adj[ida]->InsertUniqueId(idb);
      } else {
        adj[idb]->InsertUniqueId(ida);
      }
    }
  }
}

int getEdgeCount(vtkIdList** adj, int noOfVertices){
  int i, edgeCount;

  edgeCount = 0;
  for (i = 0; i < noOfVertices; ++i){
    edgeCount += adj[i]->GetNumberOfIds();
  }
  return edgeCount;
}


////////////////////////////////////////////////////////////////////

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif


