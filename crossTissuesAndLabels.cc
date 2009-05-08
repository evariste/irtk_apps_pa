#include <irtkImage.h>

char *tissues_name = NULL, *structures_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: crossTissuesAndLabels [tissueLabels] [structureLabels] [output] [offset]" << endl;
  cerr << "Images assumed in register.  Tissue labels assumed to be 1, 2 and 3 for" << endl;
  cerr << "wm, gm and csf respectively.  Add a single offset to the structure label for " << endl;
  cerr << "for gm voxels, add 2 offsets for csf voxels, leave wm voxels alone." << endl;
  cerr << "Crossed structures written to output image." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z, voxels, offset;
  irtkGreyImage tissues, structures, output;
  irtkGreyPixel *t_ptr, *s_ptr, *output_ptr;

  if (argc != 5){
    usage();
  }

  tissues_name  = argv[1];
  argc--;
  argv++;
  structures_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  offset = atoi(argv[1]);
  argc--;
  argv++;

  tissues.Read(tissues_name);
  structures.Read(structures_name);
  output.Read(tissues_name);

  x = tissues.GetX();
  y = tissues.GetY();
  z = tissues.GetZ();

  if (x != structures.GetX() || y != structures.GetY() || z != structures.GetZ() ){
    cout << "Input images must be same size." << endl;
    usage();
  }

  t_ptr      = tissues.GetPointerToVoxels();
  s_ptr      = structures.GetPointerToVoxels();
  output_ptr = output.GetPointerToVoxels();
  voxels = tissues.GetNumberOfVoxels();

  for (i = 0; i < voxels; ++i){
    if (*t_ptr > 0 && *s_ptr > 0){
      *output_ptr = *s_ptr + (*t_ptr - 1) * offset;
    } else {
      *output_ptr = 0;
    }

    ++t_ptr;
    ++s_ptr;
    ++output_ptr;
  }

  output.Write(output_name);

}
