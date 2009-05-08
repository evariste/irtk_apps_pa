/*=========================================================================
 
  Library   : project/applications
  Module    : $RCSfile: padding.cc,v $
  Authors   : (C)opyright Daniel Rueckert and Julia Schnabel 1994-2000++
              See COPYRIGHT statement in top level directory.
  Purpose   :
  Date      : $Date: 2005/02/14 18:11:50 $
  Version   : $Revision: 1.4 $
  Changes   : $Locker:  $
	      $Log: padding.cc,v $
	      Revision 1.4  2005/02/14 18:11:50  dr
	      Added usage message
	
	      Revision 1.3  2004/01/09 15:08:53  dr
	      *** empty log message ***
	
	      Revision 1.2  2003/04/17 14:10:51  dr
	      Merged branch
	
	      Revision 1.1.1.1.2.2  2002/06/11 10:02:48  dr
	      Added compatibility for MS Visual Studio
	
	      Revision 1.1.1.1.2.1  2002/05/13 13:45:59  dr
	      Minor changes
	
	      Revision 1.1.1.1  2000/10/29 15:12:17  dr
	      Imported sources
	

=========================================================================*/

#include <irtkImage.h>

void usage()
{
  cerr << "Usage: identifyregion [bigImage] [smallImage]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k;
  int i2, j2, k2;
  double x, y, z;
  int rx1, rx2, ry1, ry2, rz1, rz2;

  if (argc != 3){
    usage();
  }

  irtkGreyImage bigimage(argv[1]);
  irtkGreyImage smallimage(argv[2]);

  i = 0;
  i2 = -1;
  while (i2 < 0){
    x = i;
    bigimage.ImageToWorld(x, y, z);
    smallimage.WorldToImage(x, y, z);
    i2 = (int) round(x);
    i++;
  }
  rx1 = i - 1;

  i = bigimage.GetX();
  i2 = smallimage.GetX();
  while (i2 >= smallimage.GetX()){
    x = i;
    bigimage.ImageToWorld(x, y, z);
    smallimage.WorldToImage(x, y, z);
    i2 = (int) round(x);
    i--;
  }
  rx2 = i + 2;

  j = 0;
  j2 = -1;
  while (j2 < 0){
    y = j;
    bigimage.ImageToWorld(x, y, z);
    smallimage.WorldToImage(x, y, z);
    j2 = (int) round(y);
    j++;
  }
  ry1 = j - 1;

  j = bigimage.GetY();
  j2 = smallimage.GetY();
  while (j2 >= smallimage.GetY()){
    y = j;
    bigimage.ImageToWorld(x, y, z);
    smallimage.WorldToImage(x, y, z);
    j2 = (int) round(y);
    j--;
  }
  ry2 = j + 2;

  k = 0;
  k2 = -1;
  while (k2 < 0){
    z = k;
    bigimage.ImageToWorld(x, y, z);
    smallimage.WorldToImage(x, y, z);
    k2 = (int) round(z);
    k++;
  }
  rz1 = k -1;

  k = bigimage.GetZ();
  k2 = smallimage.GetZ();
  while (k2 >= smallimage.GetZ()){
    z = k;
    bigimage.ImageToWorld(x, y, z);
    smallimage.WorldToImage(x, y, z);
    k2 = (int) round(z);
    k--;
  }
  rz2 = k + 2;

  if (rx1 > 0)
    cout << " -Rx1 " << rx1;

  if (rx2 < bigimage.GetX())
    cout << " -Rx2 " << rx2;

  if (ry1 > 0)
    cout << " -Ry1 " << ry1;

  if (ry2 < bigimage.GetY())
    cout << " -Ry2 " << ry2;

  if (rz1 > 0)
    cout << " -Rz1 " << rz1;

  if (rz2 < bigimage.GetZ())
    cout << " -Rz2 " << rz2;

  cout  << endl;

  return 0;
}
