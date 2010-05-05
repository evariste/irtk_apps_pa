/*=========================================================================
 
  Library   : packages
  Module    : $RCSfile: dof2movie.cc,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2007-
  Purpose   :
  Date      : $Date: 2007/05/22 10:37:42 $
  Version   : $Revision: 1.1 $
  Changes   : $Locker:  $
              $Log: dof2movie.cc,v $
              Revision 1.1  2007/05/22 10:37:42  dr
              Imported sources


=========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>

// Default filenames
char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: dof2movie [dof] [n] [movie_dof_basename]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, n;
  char buffer[256];

  if (argc < 4)
    usage();
  
  // Input transformation
  dofin_name = argv[1];
  argv++;
  argc--;
  irtkTransformation *transformation = irtkTransformation::New(dofin_name);
  
  // Number of frames to generate
  n = atoi(argv[1]);
  argv++;
  argc--;
  
  // Output transformation
  dofout_name = argv[1];
  argv++;
  argc--;
  
  // Check transformation type
  irtkMultiLevelFreeFormTransformation *mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transformation);
  if (mffd != NULL){
    irtkFreeFormTransformation *ffd = mffd->GetLocalTransformation(mffd->NumberOfLevels()-1);
    
    // Store original transformation parameters
    double *params = new double[ffd->NumberOfDOFs()];
    for (j = 0; j < ffd->NumberOfDOFs(); j++){
      params[j] = ffd->Get(j);
    }
    
    // Generate sequence
    for (i = 0; i <= n; i++){
      cout << "Generating dof for t = " << (double)i / (double)n << endl;
      for (j = 0; j < ffd->NumberOfDOFs(); j++){
        ffd->Put(j, params[j] * ((double)i / (double)n));
      }
      sprintf(buffer, "%s_%.5d.dof.gz", dofout_name, i);
      cout << "Saving dof file to = " << buffer << endl;
      transformation->Write(buffer);
    }
    delete []params;
    
  } else {
    
    // Check transformation type
    irtkAffineTransformation *affine = dynamic_cast<irtkAffineTransformation *>(transformation);
    if (affine != NULL){
      
      // Store original transformation parameters
      double *params = new double[transformation->NumberOfDOFs()];
      for (j = 0; j < transformation->NumberOfDOFs(); j++){
        params[j] = transformation->Get(j);
      }
    
      // Generate sequence
      for (i = 0; i <= n; i++){
        cout << "Generating dof for t = " << (double)i / (double)n << endl;
        for (j = 0; j < transformation->NumberOfDOFs(); j++){

          if (j < 6)
            continue;

          if ((j >= 6) && (j < 9)){
            transformation->Put(j, (params[j] - 100.0) * ((double)i / (double)n) + 100.0);
          } else {
            transformation->Put(j, params[j] * ((double)i / (double)n));      
          }     
        }
        sprintf(buffer, "%s_%.5d.dof.gz", dofout_name, i);
        cout << "Saving dof file to = " << buffer << endl;
        transformation->Write(buffer);
      }
      delete []params;
    } else {
      
      // Check transformation type
      irtkRigidTransformation *rigid = dynamic_cast<irtkRigidTransformation *>(transformation);
      if (rigid != NULL){
      
        // Store original transformation parameters
        double *params = new double[transformation->NumberOfDOFs()];
        for (j = 0; j < transformation->NumberOfDOFs(); j++){
          params[j] = transformation->Get(j);
        }
    
        // Generate sequence
        for (i = 0; i <= n; i++){
          cout << "Generating dof for t = " << (double)i / (double)n << endl;
          for (j = 0; j < transformation->NumberOfDOFs(); j++){
            transformation->Put(j, params[j] * ((double)i / (double)n));
          }
          sprintf(buffer, "%s_%.5d.dof.gz", dofout_name, i);
          cout << "Saving dof file to = " << buffer << endl;
          transformation->Write(buffer);
        }
        delete []params;
      } else {
        cerr << "Not yet implemented for transformation = " << transformation->NameOfClass() << endl;
      }
    }
  }
}
