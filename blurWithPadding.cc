/*=========================================================================
 
  Library   : project/applications
  Module    : $RCSfile: blur.cc,v $
  Authors   : (C)opyright Daniel Rueckert and Julia Schnabel 1994-2000++
              See COPYRIGHT statement in top level directory.
  Purpose   :
  Date      : $Date: 2003/04/17 14:10:51 $
  Version   : $Revision: 1.2 $
  Changes   : $Locker:  $
	      $Log: blur.cc,v $
	      Revision 1.2  2003/04/17 14:10:51  dr
	      Merged branch
	
	      Revision 1.1.1.1.2.1  2002/06/11 10:03:33  dr
	      Added compatibility for MS Visual Studio
	
	      Revision 1.1.1.1  2000/10/29 15:12:17  dr
	      Imported sources
	

=========================================================================*/

#include <irtkImage.h>

#include <irtkGaussianBlurring.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: blurWithPadding [in] [out] [sigma] [pad]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkRealImage input;

  if (argc < 5){
    usage();
  }

  double sigma;
  irtkRealPixel pad = -1.0 * FLT_MAX;

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  sigma = atof(argv[1]);
  argc--;
  argv++;
  pad = atof(argv[1]);
  argc--;
  argv++;

  // Read input
  input.Read(input_name);

  // Blur image
  irtkGaussianBlurringWithPadding<irtkRealPixel> gaussianBlurring(sigma, pad);
  gaussianBlurring.SetInput (&input);
  gaussianBlurring.SetOutput(&input);
  gaussianBlurring.Run();

  // Write image
  input.Write(output_name);

  return 0;
}
