#if (defined HAS_VTK)


#include <irtkImage.h>

#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataConnectivityFilter.h>


char *input_name = NULL;

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
}

void copyFaces(vtkCellArray *src, vtkCellArray *dest){
  int i, j, noOfCells;
  vtkIdType npts = 0;
  vtkIdType *ptIds;
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
  vtkIdType *pts;

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

void printStuff(vtkPoints *points, vtkIdList** adj){

  int i, j, noOfPts;
  noOfPts = points->GetNumberOfPoints();
  double a[3];

  for (i = 0; i < noOfPts; ++i){
    points->GetPoint(i, a);
    cerr << i << " : " << a[0] << " " << a[1] << " " << a[2] << endl;
    for(j = 0; j < adj[i]->GetNumberOfIds(); ++j){
      cerr << " " << i << " - " << adj[i]->GetId(j) << endl;
    }
  }
  cerr << "Edge count: " << getEdgeCount(adj, noOfPts) << endl;

}
void usage()
{
  cerr << " polydatasphere [output] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -levels n  : Number of subdivsions, default 1." << endl;
  cerr << " -solid n   : Number to indicate which solid to start with:" << endl;
  cerr << "              0 Tetrahedron" << endl;
  cerr << "              1 Cube" << endl;
  cerr << "              2 Octahedron" << endl;
  cerr << "              3 Icosahedron (default)" << endl;
  cerr << "              4 Dodecahedron" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2){
    usage();
  }

  int ok, i;
  int noOfVerts, noOfEdges, noOfFaces;
  int noOfComponents;
  int eulerChar;
  double genus;
  int noOfArrays;

  input_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
//    if ((ok == False) && (strcmp(argv[1], "-option") == 0)){
//      argc--;
//      argv++;
//do stuff
//      ok = True;
//    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Modified();
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();
  input->Update();

  noOfVerts = input->GetNumberOfPoints();

  vtkCellArray *faces = vtkCellArray::New();
  copyFaces(input->GetPolys(), faces);

  noOfFaces = faces->GetNumberOfCells();

  vtkIdList** adj = NULL;
  createIdListArray(adj, noOfVerts);

  getAdjacency(faces, adj);
  noOfEdges = getEdgeCount(adj, noOfVerts);

  vtkPolyDataConnectivityFilter *connFilter = vtkPolyDataConnectivityFilter::New();
  connFilter->SetExtractionModeToAllRegions();
  connFilter->SetInput(input);
  connFilter->Update();

  noOfComponents = connFilter->GetNumberOfExtractedRegions();

  eulerChar = noOfVerts - noOfEdges + noOfFaces;

  genus = 0.5 * (2 * noOfComponents - eulerChar);

  cout << "V " << noOfVerts << endl;
  cout << "E " << noOfEdges << endl;
  cout << "F " << noOfFaces << endl;
  cout << "components " << noOfComponents << endl;

  cout << "Euler characteristic : V - E + F = 2C - 2g : ";
  cout << eulerChar << endl;

  cout << "Genus : " << genus << endl;

  deleteIdListArray(adj, noOfVerts);


  noOfArrays = input->GetPointData()->GetNumberOfArrays();
  if (noOfArrays > 0){
  	cout << "Arrays :" << endl;
  }
  for (i = 0; i < noOfArrays; i++){
    vtkFloatArray *scalars;

    scalars = (vtkFloatArray*) input->GetPointData()->GetArray(i);

    printf(" %2d %-20s (dim: %2d)\n", i, scalars->GetName(), scalars->GetNumberOfComponents());

  }

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

