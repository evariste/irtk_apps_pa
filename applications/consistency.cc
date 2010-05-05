/*=========================================================================
 
  Library   : packages
  Module    : $RCSfile: volumechange.cc,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2000-2004
  Purpose   :
  Date      : $Date: 2004/12/17 16:13:36 $
  Version   : $Revision: 1.1 $
  Changes   : $Locker:  $
	      $Log: volumechange.cc,v $
	      Revision 1.1  2004/12/17 16:13:36  dr
	      Imported sources
	

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

// Default filenames
char *input_name = NULL, *output_name, *dof_name1  = NULL, *dof_name2  = NULL;

void usage()
{
  cerr << "Usage: consistency [mask] [mask-value] [dof1] [dof2]\n" << endl;
  cerr << "      dof1           dof2      " << endl;
  cerr << "mask -----> image Y -----> mask" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, mask_value, xdim, ydim, zdim;
  double x1, y1, z1, x2, y2, z2, sum, sumsq, dispsq, mean;

  // Check command line
  if (argc != 5){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  mask_value = atoi(argv[1]);
  argc--;
  argv++;
  dof_name1 = argv[1];
  argc--;
  argv++;
  dof_name2 = argv[1];
  argc--;
  argv++;
  
  // Read image
  //  cout << "Reading image ... "; cout.flush();
  irtkGreyImage *image = new irtkGreyImage(input_name);
  //  cout << "done" << endl;

  // Read transformation
  irtkTransformation *transformation1;
  transformation1->irtkTransformation::Read(dof_name1);
  //  transformation1->Invert();
  irtkTransformation *transformation2;
  transformation2->irtkTransformation::Read(dof_name2);
  //  transformation2->Invert();

  n     = 0;
  sum   = 0;
  sumsq = 0;

  xdim  = image->GetX();
  ydim  = image->GetY();
  zdim  = image->GetZ();

  for (k = 0; k < zdim; k++){ 
    for (j = 0; j < ydim; j++){ 
      for (i = 0; i < xdim; i++){ 
	if (image->Get(i, j, k) == mask_value){
	  x1 = x2 = i;
	  y1 = y2 = j;
	  z1 = z2 = k;

          image->ImageToWorld(x1, y1, z1);
	  image->ImageToWorld(x2, y2, z2);

          transformation2->Transform(x2, y2, z2);
          transformation1->Transform(x2, y2, z2);

          x1 -= x2;
          y1 -= y2;
          z1 -= z2;
          dispsq = x1*x1 + y1*y1 + z1*z1;

          sum   += sqrt(dispsq);
          sumsq += dispsq;
	  n++;
	}
      }
    }
  }

  mean = sum / (double) n;
  cout <<  mean  << "," << (sumsq / (double) n) - (mean * mean) << endl;

//   cout << "mean displacement: " << mean  << endl;
//   cout << ", variance " << (sumsq / (double) n) - (mean * mean) << endl;

}
