////////////////////////////////////////////////////////////////////
// Track the distance travelled by points during Laplacian smoothing of a mesh.

#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>

char *input_name = NULL;
char *output_name = NULL;

void usage()
{
  cerr << "" << endl;
  cerr << " Usage: polydatasmooth [input] [output] [iterations] [relaxationFactor]" << endl; // <options>" << endl;
  cerr << "" << endl;
  cerr << " Area weighted Laplacian relaxation of a surface." << endl;
  cerr << "" << endl;
  cerr << " Iteratively move mesh nodes towards the centroid of the adjacent nodes." << endl;
  cerr << " The relaxation factor (RF) determines how far each node moves towards the " << endl;
  cerr << " local centroid (0<= RF <= 1).  The new position of a node is a weighted " << endl;
  cerr << " combination of its previous position and the neighbours' centroid.  " << endl;
  cerr << " RF = 1 moves a node all the way to the centroid while RF = 0 keeps a " << endl;
  cerr << " node at its previous position." << endl;
  //cerr << "Options: " << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, ok;
  int noOfIterations;
  double relaxationFactor;
  double* pts;
  double x[3];

  double currPos[3];
  unsigned short noOfCells = 0;
  vtkIdType* cells = NULL;
  vtkTriangle* triangle = NULL;
  double totalArea = 0;
  double update[3];
  vtkIdList* ptIds = NULL;
  double v1[3], v2[3], v3[3], centre[3];
  double triangleArea = 0;
  double dx, dy, dz;
  double dist, val;
  double *normal;

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
  relaxationFactor = atof(argv[1]);
  argv++;
  argc--;

  cerr << "Input        : " << input_name << endl;
  cerr << "Output       : " << output_name << endl;
  cerr << "Iterations   : " << noOfIterations << endl;
  cerr << "Relax factor : " << relaxationFactor << endl;

  // Read the polydata file
  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);

  vtkPolyData* inputRead = reader->GetOutput();
  inputRead->Update();

  // Normals need to be recalculated before using.
  vtkPolyDataNormals *normalsFilter2 = vtkPolyDataNormals::New();
  normalsFilter2->SplittingOff();
  normalsFilter2->SetInput(inputRead);
  normalsFilter2->Modified();
  normalsFilter2->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = normalsFilter2->GetOutput();

  input->Update();
  input->BuildCells();
  input->BuildLinks();

  int noOfPoints = input->GetNumberOfPoints();

  vtkFloatArray *normals = vtkFloatArray::New();
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(noOfPoints);
  normals = (vtkFloatArray*) input->GetPointData()->GetNormals();

  // the points array
  pts = new double[3*noOfPoints];

  vtkFloatArray *dists = vtkFloatArray::New();
  dists->SetNumberOfComponents(1);
  dists->SetNumberOfTuples(noOfPoints);
  for (j = 0; j < noOfPoints; ++j){
    dists->SetTuple1(j, 0.0);
  }

  // perform the smoothing
  vtkPoints* pts_original = input->GetPoints();

  // Iteration loop.
  cerr << "Iterating";
  for (i = 0; i < noOfIterations; ++i){
    cerr << ".";
    cerr.flush();

    // Loop over surface.
    for (j = 0; j < noOfPoints; ++j){

      // Initialisation for current point.
      totalArea = 0;
      cells = NULL;

      update[0] = 0;
      update[1] = 0;
      update[2] = 0;

      // Store the current position of the node.
      input->GetPoint(j, currPos);

      // What cells does this node adjoin?
      input->GetPointCells(j, noOfCells, cells);

      if ( cells == NULL )
        continue;

      for (k = 0; k < noOfCells; ++k){
        triangle = vtkTriangle::SafeDownCast(input->GetCell(cells[k]));

        if ( triangle != NULL ){
          ptIds = triangle->GetPointIds();

          input->GetPoint(ptIds->GetId(0), v1);
          input->GetPoint(ptIds->GetId(1), v2);
          input->GetPoint(ptIds->GetId(2), v3);

          triangleArea = vtkTriangle::TriangleArea(v1, v2, v3);
          vtkTriangle::TriangleCenter(v1, v2, v3, centre);

          totalArea += triangleArea;

          update[0] += triangleArea * centre[0];
          update[1] += triangleArea * centre[1];
          update[2] += triangleArea * centre[2];
        }
      }

      update[0] /= totalArea;
      update[1] /= totalArea;
      update[2] /= totalArea;

      dx = relaxationFactor * (update[0] - currPos[0]);
      dy = relaxationFactor * (update[1] - currPos[1]);
      dz = relaxationFactor * (update[2] - currPos[2]);

      dist = sqrt(dx*dx + dy*dy + dz*dz);

      normal = normals->GetTuple3(j);
      val = normal[0]*dx + normal[1]*dy + normal[2]*dz;
      if (val < 0){
        dist = -1.0 * dist;
      }

      val = dists->GetTuple1(j);
      dists->SetTuple1(j, val + dist);

      pts[j*3]   = currPos[0] + dx;
      pts[j*3+1] = currPos[1] + dy;
      pts[j*3+2] = currPos[2] + dz;
    }

    for (j = 0; j < noOfPoints; ++j){
      pts_original->SetPoint(j, pts + j*3);
    }

    input->SetPoints(pts_original);
    input->Update();
  }

  dists->SetName("smoothingDists");
  input->GetPointData()->AddArray(dists);
//  input->GetPointData()->SetActiveScalars("smoothingDists");

  // Normals need to be recalculated before saving.
  cerr << endl << "Recalculating normals";
  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();
  normalsFilter->SplittingOff();
  normalsFilter->SetInput(input);
  normalsFilter->Modified();
  normalsFilter->Update();


  // save as a vtk file
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(normalsFilter->GetOutput());
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Update();
  writer->Write();

  reader->Delete();
  writer->Delete();
  delete [] pts;
  return 0;
}




#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif



/*

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
  cerr << "Usage: tmp [input] [output] <options>\n" << endl;
  cerr << "" << endl;
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

*/

/*

/////////////////////////////////////////////////////////////
// Testing temp registration class based on fluid type approach
// with ffds.


#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>

char *target_name = NULL;
char *source_name = NULL;
char *dofin_name = NULL;
char *dofout_name = NULL;
char *parameter_name = NULL;
char *mask_name = NULL;

void usage()
{
}

int main(int argc, char **argv)
{
//   int iterations = 1;
//   int ok;
//   double dx, dy, dz;
//   int paddingValue;

//   // Parse source and target images
//   target_name = argv[1];
//   argc--;
//   argv++;
//   source_name = argv[1];
//   argc--;
//   argv++;

//   cout << "Reading target and source ... "; cout.flush();
//   irtkGreyImage target(target_name);
//   irtkGreyImage source(source_name);
//   cout << "done" << endl;

//   dx = 20;
//   dy = 20;
//   dz = 20;

//   while (argc > 1){
//     ok = False;
//     if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)){
//       argc--;
//       argv++;
//       dofin_name = argv[1];
//       argc--;
//       argv++;
//       ok = True;
//     }
//     if ((ok == False) && (strcmp(argv[1], "-dofout") == 0)){
//       argc--;
//       argv++;
//       dofout_name = argv[1];
//       argc--;
//       argv++;
//       ok = True;
//     }
//     if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
//       argc--;
//       argv++;
//       paddingValue = atoi(argv[1]);
//       argc--;
//       argv++;
//       ok = True;
//     }
//     if ((ok == False) && (strcmp(argv[1], "-ds") == 0)){
//       argc--;
//       argv++;
//       dx = atof(argv[1]);
//       dy = atof(argv[1]);
//       dz = atof(argv[1]);
//       argc--;
//       argv++;
//       ok = True;
//     }
//     if ((ok == False) && (strcmp(argv[1], "-parameter") == 0)){
//       argc--;
//       argv++;
//       ok = True;
//       parameter_name = argv[1];
//       argc--;
//       argv++;
//     }
//     if ((ok == False) && (strcmp(argv[1], "-mask") == 0)){
//       argc--;
//       argv++;
//       mask_name = argv[1];
//       argc--;
//       argv++;
//       ok = True;
//     }
//     if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)){
//       argc--;
//       argv++;
//       iterations = atoi(argv[1]);
//       argc--;
//       argv++;
//       ok = True;
//     }
//     if (ok == False){
//       cerr << "Can not parse argument " << argv[1] << endl;
//       usage();
//     }
//   }

//   // Create transformation
//   irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;

//   // Read transformation
//   if (dofin_name != NULL){
//     irtkTransformation *transform = irtkTransformation::New(dofin_name);
//     if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0){
//       mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
//     } else {
//       if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0){
// 	mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
//       } else {
// 	if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0){
// 	  mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
// 	} else {
// 	  cerr << "Input transformation is not of type rigid, affine " << endl;
//  	  cerr << "or multi-level free form deformation" << endl;
// 	  exit(1);
// 	}
//       }
//     }
//     delete transform;
//   } else {
//     mffd = new irtkMultiLevelFreeFormTransformation;
//   }

//   // Create ffd
//   irtkBSplineFreeFormTransformation *affd = new irtkBSplineFreeFormTransformation(target, dx, dy, dz);

//   // Add ffd
//   mffd->PushLocalTransformation(affd);



//   irtkTempRegistration *registration = new irtkTempRegistration;

//   // Set input and output for the registration filter
//   registration->SetInput(&target, &source);
//   registration->SetOutput(mffd);

//   registration->SetIterations(iterations);

//   // Read default parameter
//   registration->irtkImageRegistration::Read(parameter_name);

//   // Run registration filter
//   registration->Run();

//   mffd->irtkTransformation::Write(dofout_name);

}


*/


/*

/////////////////////////////////////////////////////////////
/// Can we find the rigid components of the Jacobian tensor?

#include <irtkImage.h>

#include <irtkTransformation.h>

// Default filenames
char *input_name = NULL, *output_name, *dof_name  = NULL;

typedef enum { LocalJacobian, GlobalJacobian, RelativeJacobian, TotalJacobian } JacobianMode;

void usage()
{
  cerr << "Usage: jacobian [input] [output] [ffd]\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-total>                 Total jacobian (default)" << endl;
  cerr << "<-local>                 Local jacobian only " << endl;
  cerr << "<-global>                Global jacobian only " << endl;
  cerr << "<-relative>              Local jacobian divided by global Jacobian" << endl;
  cerr << "<-padding value>         Padding value" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, m, ok, fluid, padding;
  double x, y, z;
  irtkMatrix jac33(3, 3);
  irtkMatrix jac44(4, 4);
  irtkMatrix tmpmat;

  jac44.Ident();

  irtkMatrix L, M, N;
  irtkVector w;

  irtkAffineTransformation *trans = new irtkAffineTransformation;

  // Check command line
  if (argc < 4){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Initialize padding value
  padding = MIN_GREY;

  // Initialize jacobian mode
  fluid = False;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkGreyImage *image = new irtkGreyImage(input_name);
  cout << "done" << endl;

  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd;
  if (fluid == True){
    mffd = new irtkFluidFreeFormTransformation;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }
  mffd->irtkTransformation::Read(dof_name);

  for (k = 0; k < image->GetZ(); k++){
    for (j = 0; j < image->GetY(); j++){
      for (i = 0; i < image->GetX(); i++){
	if (image->Get(i, j, k) > padding){

	  x = i;
	  y = j;
	  z = k;
	  image->ImageToWorld(x, y, z);
          mffd->Jacobian(x, y, z, jac33);

          for (m = 0; m < 3; ++m){
            for (n = 0; n < 3; ++n){
              jac44(m, n) = jac33(m, n);
            }
          }

          trans->PutMatrix(jac44);
          trans->Print();
          cout << endl;

          // Studholme approach:

          tmpmat = jac33;
          tmpmat.Transpose();
          tmpmat = jac33 * tmpmat;
          tmpmat = sqrtm(tmpmat);
          tmpmat.Invert();
          tmpmat = tmpmat * jac33; // ( J * J^T )^(-1/2) * J

          for (m = 0; m < 3; ++m){
            for (n = 0; n < 3; ++n){
              jac44(m, n) = tmpmat(m, n);
            }
          }
          trans->PutMatrix(jac44);
          trans->Print();
          cout << endl;


          tmpmat = jac33;
          tmpmat.Transpose();
          tmpmat = tmpmat * jac33;
          tmpmat = sqrtm(tmpmat);
          tmpmat.Invert();
          tmpmat = jac33 * tmpmat; // J * (J^T * J)^(-1/2)

          for (m = 0; m < 3; ++m){
            for (n = 0; n < 3; ++n){
              jac44(m, n) = tmpmat(m, n);
            }
          }
          trans->PutMatrix(jac44);
          trans->Print();
          cout << endl;



          tmpmat.Invert();
          tmpmat = tmpmat * jac33;
          for (m = 0; m < 3; ++m){
            for (n = 0; n < 3; ++n){
              jac44(m, n) = tmpmat(m, n);
            }
          }

          trans->PutMatrix(jac44);
          trans->Print();

          cout << endl;

          // Guimond approach:

          L.Initialize(3, 3);
          M.Initialize(3, 3);
          N.Initialize(3, 3);
          w.Initialize(3);

          jac33.SVD(L, w, N);
          // jac33 = L * M * N^T
          M(0, 0) = w(0);
          M(1, 1) = w(1);
          M(2, 2) = w(2);

          N.Transpose();
          tmpmat = L * N;

          for (m = 0; m < 3; ++m){
            for (n = 0; n < 3; ++n){
              jac44(m, n) = tmpmat(m, n);
            }
          }

          trans->PutMatrix(jac44);
          trans->Print();
          cout << endl;

          tmpmat.Invert();
          tmpmat = tmpmat * jac33;
          for (m = 0; m < 3; ++m){
            for (n = 0; n < 3; ++n){
              jac44(m, n) = tmpmat(m, n);
            }
          }

          trans->PutMatrix(jac44);
          trans->Print();

          cout << endl;
          tmpmat.Print();

          cout << "===============================================================================" << endl;


        }
      }
    }
  }

}

*/

