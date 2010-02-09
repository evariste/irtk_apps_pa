/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: defmap.cc 117 2009-12-11 14:40:13Z pa100 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-12-11 14:40:13 +0000 (Fri, 11 Dec 2009) $
  Version   : $Revision: 117 $
  Changes   : $Author: pa100 $

=========================================================================*/

#ifdef HAS_VTK

#include <irtkRegistration.h>
#include <nr.h>

// Default filenames
char *target_name = NULL, *output_name = NULL, *mask_name = NULL;
char **dof_name  = NULL;

#define MAX_DOFS 10

#define MAX_PTS_PAREG 10000

void usage()
{
  cerr << " Usage: defmapstats [target] <options>" << endl;
  cerr << " " << endl;
  cerr << " Stats of distances that voxels in the target move under the effect" << endl;
  cerr << " of one or more transformations." << endl;
  cerr << " " << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " -dofin file      Transformation. Can be repeated to give multiple dofs" << endl;
  cerr << "                  (i.e. \"-dofin file1 -dofin file2 ...\") which are applied" << endl;
  cerr << "                  in order given (Max 10)." << endl;
  cerr << " -Tp  value       Padding value in target." << endl;
  cerr << " -mask file       Image mask to define a region of interest." << endl;
  cerr << " -removeGlobal    Estimate a global affine transformation based on a sampled" << endl;
  cerr << "                  subset of voxel locations and their transformed coordinates." << endl;
  cerr << "                  Remove the effect of this global transformation before " << endl;
  cerr << "                  estimating the distances." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation **transformation = NULL;
  int ok, x, y, z, regressAffine, squaredDistance;
  double Tp, val;
  int i, noOfDofs;
  int noOfPoints, ptID;
  double sumD, sumD2;
  

  // Check command line
  if (argc < 3){
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  
  // Read target image
  irtkRealImage target(target_name);

  Tp = -1.0 * FLT_MAX;
  
  // Fix number of dofs
  noOfDofs = 0;

  dof_name = new char*[MAX_DOFS];
  
  regressAffine = False;
  squaredDistance = False;

  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)){
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      noOfDofs++;
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)){
      argc--;
      argv++;
      Tp = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-removeGlobal") == 0)){
      argc--;
      argv++;
      regressAffine = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-square") == 0)){
      argc--;
      argv++;
      squaredDistance = True;
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
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  if (noOfDofs == 0){
    cerr << "No transformations specified: Using a single identity transformation!" << endl;
    noOfDofs = 1;
    transformation = new irtkTransformation*[noOfDofs];
    transformation[0] = new irtkRigidTransformation;
  } else {
    // Read dof(s)
    transformation = new irtkTransformation*[noOfDofs];
    for (i = 0; i < noOfDofs; i++){
      transformation[i] = irtkTransformation::New(dof_name[i]);
    }
  }
  
  // Setting up mask
  cout.flush();
  
  irtkGreyImage mask;
  
  if (mask_name == NULL){
    mask.Read(target_name);
    irtkGreyPixel *ptr2mask = mask.GetPointerToVoxels();
    irtkRealPixel *ptr2tgt  = target.GetPointerToVoxels();
    noOfPoints = target.GetNumberOfVoxels();
    
    for (i = 0; i < noOfPoints; i++){
      if (*ptr2tgt > Tp)
        *ptr2mask = 1;
      else
        *ptr2mask = 0;
      
      ++ptr2tgt;
      ++ptr2mask;
    }
    
  } else {
    mask.Read(mask_name);
  }

  // Make an identity global transformation.
  irtkAffineTransformation *trAffine = new irtkAffineTransformation;
  
  
  if (regressAffine == True){
    // Estimate the global affine transformation.
    
    irtkPointSet targetPts;
    irtkPointSet sourcePts;
    
    // Collect point data.

    noOfPoints = target.GetNumberOfVoxels();

    int incr;
    incr = 1;
    while ((noOfPoints / incr) > MAX_PTS_PAREG){
      incr++;
    }
    
    // Subsample uniformly by increments of 'incr' 
    noOfPoints = 0;
    ptID = -1;

    // Loop over all voxels.
    for (z = 0; z < target.GetZ(); z++){
      for (y = 0; y < target.GetY(); y++){
        for (x = 0; x < target.GetX(); x++){

          ptID++;
          // Should we sample it or not?
          if ((ptID % incr) != 0)
            continue;

          // Get two copies of current image coordinates.
          irtkPoint p(x, y, z);
          irtkPoint q(x, y, z);

          // Transform points into target world coordinates.
          target.ImageToWorld(p);
          target.ImageToWorld(q);

          targetPts.Add(p);

          // Transform one point to source coordinates.
          for (i = 0; i < noOfDofs; i++){
            transformation[i]->Transform(q);
          }

          sourcePts.Add(q);

          noOfPoints++;
        }
      }
    }

    
    // Estimate global affine component
    
    irtkPointAffineRegistration *pareg = new irtkPointAffineRegistration;
    // Set input and output for the registration filter
    irtkPointSet tmp1 = targetPts;
    irtkPointSet tmp2 = sourcePts;
    pareg->SetInput(&tmp1, &tmp2);
    pareg->SetOutput(trAffine);

    // Run registration filter
    pareg->Run();

    
  }
  
  // Finally can calculate distances.
  noOfPoints = 0;
  for (z = 0; z < target.GetZ(); z++){
    for (y = 0; y < target.GetY(); y++){
      for (x = 0; x < target.GetX(); x++){

        if (mask(x,y,z) > 0){
          noOfPoints++;
        }
      }
    }
  }

  float *vals, *vals2;
  vals  = new float[1 + noOfPoints];
  vals2 = new float[1 + noOfPoints];
  
  noOfPoints = 0;
  for (z = 0; z < target.GetZ(); z++){
    for (y = 0; y < target.GetY(); y++){
      for (x = 0; x < target.GetX(); x++){

        if (mask(x,y,z) > 0){

          irtkPoint p(x, y, z);
          irtkPoint q(x, y, z);

          // Transform points into target world coordinates
          target.ImageToWorld(p);
          target.ImageToWorld(q);
          
          // Apply global Affine to one copy of the target 
          // points (this is the identity if no regression was done).
          trAffine->irtkTransformation::Transform(p);
          
          // Transform the other copy by the input dofs.
          for (i = 0; i < noOfDofs; i++){
            transformation[i]->Transform(q);
          }
          
          // Calculate distance.
          val = p.Distance(q);
          sumD += val;
          sumD2 += val * val;

          vals[1 + noOfPoints]  = val;
          vals2[1 + noOfPoints] = val * val;
          noOfPoints++;
        }

      }
    }
  }


  sort(noOfPoints, vals);
  sort(noOfPoints, vals2);

  i = 1 + (int) round(0.5 * (noOfPoints - 1));

  cout << endl;
  cout << sumD << " " << sumD2 << " " ;
  cout << noOfPoints << " " << sumD / (double (noOfPoints));
  cout << " " << sumD2 / (double (noOfPoints)) << " ";
  cout << vals[i] << " " << vals2[i] << endl;
  
//      index = 1 + (int) round( (double) percentile[i] * (count - 1) / 100.0);
    
    
  
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
