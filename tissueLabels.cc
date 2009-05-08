#include <irtkImage.h>

char *wm_name = NULL, *gm_name = NULL, *csf_name = NULL, *output_name = NULL;

#define MAX(a,b){((a)>(b)) ? (a) : (b);}


void usage()
{
  cerr << "Usage: tissueLabels [wm] [gm] [csf] [out]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z, voxels, max, tieCount = 0;
  irtkGreyImage wm, gm, csf, output;
  irtkGreyPixel *wm_ptr, *gm_ptr, *csf_ptr, *output_ptr;

  if (argc != 5){
    usage();
  }

  wm_name  = argv[1];
  argc--;
  argv++;
  gm_name  = argv[1];
  argc--;
  argv++;
  csf_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  wm.Read(wm_name);
  gm.Read(gm_name);
  csf.Read(csf_name);
  output.Read(wm_name);

  x = wm.GetX();
  y = wm.GetY();
  z = wm.GetZ();

  if (x != gm.GetX() || x != csf.GetX() || y != gm.GetY() || y != csf.GetY() || z != gm.GetZ() || z != csf.GetZ()){
    cout << "Input images must be same size." << endl;
    usage();
  }

  wm_ptr     = wm.GetPointerToVoxels();
  gm_ptr     = gm.GetPointerToVoxels();
  csf_ptr    = csf.GetPointerToVoxels();
  output_ptr = output.GetPointerToVoxels();
  voxels = wm.GetNumberOfVoxels();

  for (i = 0; i < voxels; ++i){
    if (*wm_ptr > 0 || *gm_ptr > 0 || *csf_ptr > 0){

      max = *wm_ptr;
      if (max < *gm_ptr)
        max = *gm_ptr;
      if (max < *csf_ptr)
        max = *csf_ptr;

      if (*wm_ptr == *gm_ptr && *wm_ptr == *csf_ptr){
        ++tieCount; // all 3 equal
      } else if(max == *wm_ptr && max == *gm_ptr){
        ++tieCount; // two tied for max
      } else if(max == *wm_ptr && max == *csf_ptr){
        ++tieCount; 
      } else if(max == *gm_ptr && max == *csf_ptr){
        ++tieCount;
      }

      if (max == *wm_ptr){
        *output_ptr = 1;
      } else if (max == *gm_ptr){
        *output_ptr = 2;
      } else if (max == *csf_ptr){
        *output_ptr = 3;
      }

    } else { // all tissue maps are zero.
      *output_ptr = 0;
    }

    ++wm_ptr;
    ++gm_ptr;
    ++csf_ptr;
    ++output_ptr;
  }

  output.Write(output_name);

  cout << "There were " << tieCount << " ties" << endl;

}
