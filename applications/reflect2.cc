#include <irtkImage.h>
#include <irtkFileToImage.h>

char *input_name = NULL, *output_name = NULL;

char **reflect_list = NULL;

void usage()
{
  cerr << "Usage: reflect [in] [out] [reflection_1] <reflection_2> .. <reflection_n>" << endl;
  cerr << "Where the reflections are chosen from:" << endl;
  cerr << " " << endl;
  cerr << "        -x -y -z -xy -xz -yz" << endl;
  cerr << " " << endl;
  cerr << "The reflections are processed in the order given." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  irtkBaseImage *image;
  
  if (argc < 4){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  int i;
  int no = argc - 1;
  reflect_list = new char *[no];
  for (i = 0; i < no; ++i){
    reflect_list[i] = argv[1];
    argc--;
    argv++;
  }

  // Read input
  irtkFileToImage *reader = irtkFileToImage::New(input_name);
  image = reader->GetOutput();
  
  for (i = 0; i < no; ++i){
    if (strcmp("-x", reflect_list[i]) == 0){
      image->ReflectX();
    }

    if (strcmp("-y", reflect_list[i]) == 0){
     image->ReflectY();
    }

    if (strcmp("-z", reflect_list[i]) == 0){
        image->ReflectZ();
    }

    if (strcmp("-xy", reflect_list[i]) == 0 || 
        strcmp("-yx", reflect_list[i]) == 0){
      image->FlipXY(0);
    }

    if (strcmp("-xz", reflect_list[i]) == 0 ||
        strcmp("-zx", reflect_list[i]) == 0){
      image->FlipXZ(0);
    }

    if (strcmp("-yz", reflect_list[i]) == 0 ||
        strcmp("-zy", reflect_list[i]) == 0){
      image->FlipYZ(0);
    }
  }
 
  // Write region
  image->Write(output_name);

  delete [] reflect_list;

  return 0;
}
