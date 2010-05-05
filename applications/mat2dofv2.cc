/*=========================================================================
 
  Library   : packages
  Module    : $RCSfile: mat2dof.cc,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2000-2001
  Purpose   :
  Date      : $Date: 2003/04/17 14:10:51 $
  Version   : $Revision: 1.2 $
  Changes   : $Locker:  $
	      $Log: mat2dof.cc,v $
	      Revision 1.2  2003/04/17 14:10:51  dr
	      Merged branch
	
	      Revision 1.1.1.1.2.2  2003/01/03 17:21:12  dr
	      Imported sources
	

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

void usage()
{
  cerr << "Usage: mat2dof [matfile] [doffile] [-invert]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkMatrix matrix;
  irtkAffineTransformation transformation;

  if ((argc != 3) && (argc != 4)){
    usage();
  }

  // Read matrix
  matrix.Read(argv[1]);
  matrix.Print();

  if (argc == 4){
    if (strcmp(argv[3], "-invert") == 0){
      matrix.Invert();
    } else {
      usage();
    }
  }

  // Convert to rigid transformation
  transformation.PutMatrix(matrix);

  // Write transformation
  transformation.irtkTransformation::Write(argv[2]);
  transformation.Print();
}
