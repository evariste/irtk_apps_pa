#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

char *input_name  = NULL;
char *output_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << " " << endl;
  cerr << " Usage: polydatascalarsmooth [input] [output] [iterations] [kernel]" << endl; // <options>" << endl;
  cerr << "" << endl;
  cerr << " Smooth the scalar data associated with a surface." << endl;
  cerr << " Currently only smooths the active scalars." << endl;
  cerr << " " << endl;
  cerr << " Use Laplacian smoothing (local averaging) to find a smoothed estimate for" << endl;
  cerr << " the scalar data.  The average given to a point is obtained from the scalar" << endl;
  cerr << " value for the point and those of each of its edge neighbours.  The weights" << endl;
  cerr << " are one for the current point and a Gaussian type weight for the edge" << endl;
  cerr << " neighbours that diminishes with the length of the edge." << endl;
  cerr << " " << endl;
  cerr << " This is carried out iteratively.  More iterations and / or a larger kernel" << endl;
  cerr << " result in smoother output." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " <-scalarName>    Name of scalars to smooth." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  vtkIdType i;
  long j, k, n, ok;
  int noOfIterations;
  double kernel;
	int noOfPoints, noOfCells;

	double v1[3], v2[3];
  vtkIdType cellType;
  vtkIdType *ptsInCell = NULL;
  vtkIdType noOfPointsInCell;
  int u, v;

  int noOfEdgesAtPt;
  double sumDist, dist, meanDist;
  int count;
  double val, sumVals, distSq;
  double w, sumW;

  if (argc < 4){
    usage();
  }

  // Parse arguments.
  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;
  noOfIterations = atoi(argv[1]);
  argv++;
  argc--;
  kernel = atof(argv[1]);
  argv++;
  argc--;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-scalarName") == 0)){
      argc--;
      argv++;
      scalar_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }

    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      exit(1);
    }
  }

  cerr << "Input      : " << input_name << endl;
  cerr << "Output     : " << output_name << endl;
  cerr << "Iterations : " << noOfIterations << endl;
  cerr << "Sigma      : " << kernel << endl;

  if (scalar_name != NULL){
    cerr << "Scalars  : " << scalar_name << endl;
  } else {
    cerr << "Scalars  : default" << endl;
  }

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData* input = reader->GetOutput();
  input->BuildCells();
  input->BuildLinks();
  input->Update();

	noOfPoints = input->GetNumberOfPoints();
  noOfCells  = input->GetNumberOfCells();

  // Storage for edge information.  Really this is double the amount needed
  // as the edge information for each edge's endpoint are stored.
  vtkIdList **edges;

  edges = new vtkIdList *[noOfPoints];
  for (i = 0; i < noOfPoints; ++i){
    edges[i] = vtkIdList::New();
  }
  cerr << "Allocated edges " << endl;

  // Find the adjacency information.
  for (i = 0; i < noOfCells; ++i){
    input->GetCellPoints(i, noOfPointsInCell, ptsInCell);

    u = ptsInCell[0];
    v = ptsInCell[noOfPointsInCell-1];

    edges[u]->InsertUniqueId(v);
    edges[v]->InsertUniqueId(u);

    for (j = 0; j < noOfPointsInCell - 1 ; ++j){
      u = ptsInCell[j];
      v = ptsInCell[j + 1];

      edges[u]->InsertUniqueId(v);
      edges[v]->InsertUniqueId(u);
    }
  }

  // Calculate mean edge length. The following is a bit redundant as every edge will
  // appear twice in the calculation.
  sumDist = 0.0;
  count   = 0;

  for (i = 0; i < noOfPoints; ++i){
    noOfEdgesAtPt = edges[i]->GetNumberOfIds();

    input->GetPoint(i, v1);

    for (j = 0; j < noOfEdgesAtPt; ++j){
      input->GetPoint(edges[i]->GetId(j), v2);

      dist = sqrt(vtkMath::Distance2BetweenPoints(v1, v2));
      sumDist += dist;
      ++count;
    }
  }

  meanDist = sumDist / ((double) count);

  cerr << "Points      : " << noOfPoints << endl;
  cerr << "Cells       : " << noOfCells << endl;
  cerr << "Edges       : " << count / 2 << endl;
  cerr << "Mean length : " << meanDist << endl;

  // Storage for the scalars at the start of each iteration.
  vtkFloatArray *scalarsIn = vtkFloatArray::New();
  // Storage for the updated scalars.
  vtkFloatArray *scalarsOut = vtkFloatArray::New();

  // Initialise the updated scalars' size by copying the current.
  if (scalar_name != NULL){
    scalarsOut = (vtkFloatArray *) input->GetPointData()->GetScalars(scalar_name);
  } else {
    scalarsOut = (vtkFloatArray *) input->GetPointData()->GetScalars();
  }
  cerr << "Iterating"; cerr.flush();

  double denom = 2 * kernel * kernel * meanDist * meanDist;
  
  for (n = 0; n < noOfIterations; ++n){
    cerr << "."; cerr.flush();

    // Get the current scalars.
    if (scalar_name != NULL){
      scalarsIn = (vtkFloatArray *) input->GetPointData()->GetScalars(scalar_name);
    } else {
      scalarsIn = (vtkFloatArray *) input->GetPointData()->GetScalars();
    }

    for (i = 0; i < noOfPoints; ++i){
      input->GetPoint(i, v1);

      noOfEdgesAtPt = edges[i]->GetNumberOfIds();
      // Initialise by including the current point.
      sumVals = scalarsIn->GetTuple1(i);
      sumW    = 1.0;

      for (j = 0; j < noOfEdgesAtPt; ++j){
        k = edges[i]->GetId(j);
        input->GetPoint(k, v2);

        distSq = vtkMath::Distance2BetweenPoints(v1, v2);
        w = exp(-1.0 * distSq / denom);

        // Add adjoining point's contribution.
        sumVals += w * scalarsIn->GetTuple1(k);
        sumW    += w;
      }
      // Updated value for current point.
      scalarsOut->SetTuple1(i, sumVals / sumW);
    }

    if (scalar_name != NULL){
      input->GetPointData()->AddArray(scalarsOut);
    } else {
      input->GetPointData()->SetScalars(scalarsOut);
    }
      input->Update();
  }

	// Write the result.
	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInput(input);
	writer->SetFileName(output_name);
	writer->SetFileTypeToBinary();
  writer->Update();
	writer->Write();

  // Clean up.
	reader->Delete();
	writer->Delete();

  for (i = 0; i < noOfPoints; ++i){
    edges[i]->Delete();
  }
  delete [] edges;

	return 0;
}




#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

/*

  // INSERTING EDGES UNIQUELY:

  for (i = 0; i < noOfCells; ++i){

    cellType = input->GetCellType(i);
    if (cellType == VTK_POLYGON || cellType == VTK_TRIANGLE || cellType == VTK_LINE){
      input->GetCellPoints(i, noOfPointsInCell, ptsInCell);

      u = ptsInCell[0];
      v = ptsInCell[noOfPointsInCell-1];

      if (u < v){
        edges[u]->InsertUniqueId(v);
      } else if (v < u){
        edges[v]->InsertUniqueId(u);
      }

      for (j = 0; j < noOfPointsInCell - 1 ; ++j){
        u = ptsInCell[j];
        v = ptsInCell[j + 1];

        if (u < v){
          edges[u]->InsertUniqueId(v);
        } else if (v < u){
          edges[v]->InsertUniqueId(u);
        }
      }

    }
  }


*/
