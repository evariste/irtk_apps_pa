
// Firstly, read the FFD file, loop over the control points, and multiply
// them by the affine component so that we have replaced the FFD.

// Then, loop over the reference surface, and transform the point into the
// source space using just the mffd part of the transformation and replace
// the surface point. Ouput this surface.

#ifdef HAS_VTK

#include <irtkImage.h>
#include <irtkTransformation.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

char *in_surface_name = NULL, *out_surface_name = NULL, *dof_name  = NULL;

void usage(){
  cout << "Usage: stransformation_local [surface in] [surface out] <-dofin dofFile>" << endl;
  cout << "\nTransform [surface in] using the local component of given transformation." << endl;
  cout << "The local component is affine corrected first." << endl;
  cout << "Result written to [surface out]\n" << endl;
  exit(0);
}

int main(int argc, char **argv)
{
  int i, j, k;
  bool ok;
  double p[3];

  if (argc < 3){
    usage();
  }

  // Parse input and output surface names
  in_surface_name  = argv[1];
  argc--;
  argv++;
  out_surface_name = argv[1];
  argc--;
  argv++;


  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  irtkTransformation *transform;
  irtkMultiLevelFreeFormTransformation *mffd;
  irtkBSplineFreeFormTransformation *affd;

  if (dof_name != NULL){
    // Read transformation
    transform = irtkTransformation::New(dof_name);
    mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);
  } else {
    // Create identity transformation
    mffd = new irtkMultiLevelFreeFormTransformation;
    affd = new irtkBSplineFreeFormTransformation;
    mffd->PushLocalTransformation(affd);
  }

  // Lattice dimensions.
  int xdim = mffd->GetLocalTransformation(0)->GetX();
  int ydim = mffd->GetLocalTransformation(0)->GetY();
  int zdim = mffd->GetLocalTransformation(0)->GetZ();

  double corrected_disp[3];
  double disp[3];

  // Get the global transformation.
  irtkMatrix m = mffd->GetMatrix();
  // The control point displacements are corrected by multiplying by the
  // Jacobian of the global transformation (upper left 3x3 block) from
  // source image to target image.  The matrix m therefore needs to be
  // inverted (as it represents the target to source global
  // transformation).
  m.Invert();

  for (k = 0; k < zdim; k++){
    for (j = 0; j < ydim; j++){
      for (i = 0; i < xdim; i++){
        mffd->GetLocalTransformation(0)->Get(i, j, k, disp[0], disp[1], disp[2]);

        corrected_disp[0] = m(0, 0) * disp[0] + m(0, 1) * disp[1] + m(0, 2) * disp[2];
        corrected_disp[1] = m(1, 0) * disp[0] + m(1, 1) * disp[1] + m(1, 2) * disp[2];
        corrected_disp[2] = m(2, 0) * disp[0] + m(2, 1) * disp[1] + m(2, 2) * disp[2];

        mffd->GetLocalTransformation(0)->Put(i, j, k, corrected_disp[0], corrected_disp[1], corrected_disp[2]);
      }
    }
  }

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(in_surface_name);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

  for (i = 0; i < surface->GetNumberOfPoints(); i++){
    surface->GetPoints()->GetPoint(i, p);

    // Transform the point using the local FFD only.
    mffd->GetLocalTransformation(0)->Transform(p[0], p[1], p[2]);
  
    surface->GetPoints()->SetPoint(i, p);
 }

  surface->Modified();
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(out_surface_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library"
       << endl;
}
#endif
