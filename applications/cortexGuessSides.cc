/////////////////////////////////////////////////////////////////////////
// TRY AND IDENTIFY THE INNER AND OUTER SURFACES OF AN EXTRACTED CORTICAL
// SURFACE BASED ON NORMALS AND THE CENTRE OF GRAVITY.

#if (defined HAS_VTK)


#include <irtkImage.h>

#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: cortexGuessSides [input] [output] <options>\n" << endl;
  cerr << "The input is a surface extracted from a WM map." << endl;
  cerr << "The outer boundary of this surface is the inner cortical surface." << endl;
  cerr << "The inner boundary of this surface is not so interesting." << endl;

  cerr << "Try and guess which is the inner and outer surfaces based" << endl;
  cerr << "on the relative direction of the normal at a point and the" << endl;
  cerr << "displacement vector from the point to the centre of gravity" << endl;
  cerr << "of the whole surface." << endl;

  cerr << "Assign scalar values to indicate the cosine of the angle " << endl;
  cerr << "between these two vectors." << endl;
  cerr << "" << endl;
  cerr << "-threshold : Assign 0 or 1 based on the sign of the cosine." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  int threshold = false;

  if (argc < 3){
    usage();
  }

  // Parse source and target point lists
  in_name  = argv[1];
  argc--;
  argv++;
  out_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;
      threshold = true;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read the input surface.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();

  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();

  normalsFilter->SetInput(input);
  normalsFilter->SplittingOff();
  normalsFilter->Modified();
  normalsFilter->Update();

  vtkPolyData *surf = normalsFilter->GetOutput();
  surf->Update();


  double cx, cy, cz;
  int i, n;
  double pt[3];
  double pt2cent[3];
  double *normal;
  double val;

  n = surf->GetNumberOfPoints();
  cx = cy = cz = 0.0;

  for (i = 0; i < n; i++){
    surf->GetPoint (i, pt);
    cx += pt[0];
    cy += pt[1];
    cz += pt[2];
  }

  if (n > 0){
    cx /= n;
    cy /= n;
    cz /= n;
  }

  vtkFloatArray *cosines = vtkFloatArray::New();

  cosines->SetNumberOfComponents(1);
  cosines->SetNumberOfTuples(n);

  vtkFloatArray *normals = vtkFloatArray::New();
  normals = (vtkFloatArray*) surf->GetPointData()->GetNormals();


  for (i = 0; i < n; i++){
    surf->GetPoint (i, pt);
    normal = normals->GetTuple(i);

    pt2cent[0] = cx - pt[0];
    pt2cent[1] = cy - pt[1];
    pt2cent[2] = cz - pt[2];

    val = vtkMath::Dot(pt2cent, normal) / vtkMath::Norm(pt2cent) / vtkMath::Norm(normal);

    if (threshold == true)
      val = val < 0 ? 1 : 0;

    cosines->SetTuple1(i, val);
  }

  cosines->SetName("Cosines");
  surf->GetPointData()->AddArray(cosines);
  surf->GetPointData()->SetActiveScalars("Cosines");
  surf->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(out_name);
  writer->SetInput(surf);
  writer->Modified();
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

