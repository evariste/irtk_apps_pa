

#include "irtkImage.h"

#include "zeta.h"

char *input_name = NULL, *output_name = NULL;
char *mask_name = NULL;
char **ref_name;




void usage()
{
    cout << "Usage" << endl;
    exit(0);
}


int main(int argc, char **argv)
{

    int ok;

    int patchRadius = 0; // Diam = 1 + 2*R
    int nbhdRadius = 4;
    int refCount = 0;

    if (argc < 2){
        usage();
    }

    input_name  = argv[1];
    argc--;
    argv++;
    output_name = argv[1];
    argc--;
    argv++;

    irtkRealImage image;
    image.Read(input_name);


    while (argc > 1) {
      ok = false;

      if ((ok == false) && (strcmp(argv[1], "-patchRadius") == 0)) {
        argc--;
        argv++;
        patchRadius = atoi(argv[1]);
        argc--;
        argv++;

        cout << "Setting patch size to " << patchRadius << endl;
        ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-nbhdRadius") == 0)) {
        argc--;
        argv++;
        nbhdRadius = atoi(argv[1]);
        argc--;
        argv++;

        cout << "Setting neighbourhood size to " << nbhdRadius << endl;
        ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
        argc--;
        argv++;
        mask_name = argv[1];
        argc--;
        argv++;

        cout << "Using mask file: " << mask_name << endl;
        ok = true;
      }

      if ((ok == false) && (strcmp(argv[1], "-refs") == 0)) {
        argc--;
        argv++;

        int n = 1;
        while ((n < argc) && (argv[n][0] != '-')){
          n++;
        }

        refCount = n - 1;
        ref_name = new char*[refCount];

        for (int i = 0; i < refCount; i++){
          ref_name[i] = argv[1];
          argc--;
          argv++;
        }

        ok = true;
      }


      if (ok == false) {
        cerr << "Unknown option: " << argv[1] << endl;
        usage();
      }
    }


    if (refCount == 0){
      cout << "No reference images. Nothing to do." << endl;
      cout << "Exiting" << endl;
      exit(0);
    }

    irtkRealImage **refImg = new irtkRealImage*[refCount];

    cout << "Using the reference images: " << endl;
    for (int i = 0; i < refCount; i++){
      cout << "     " << ref_name[i] << endl;
      refImg[i] = new irtkRealImage(ref_name[i]);
    }





    Zeta zetaFilt;

    zetaFilt.SetTarget(&image);

    if (mask_name != NULL){
      irtkGreyImage *mask = new irtkGreyImage(mask_name);
      zetaFilt.SetMask(mask);
    }

    zetaFilt.SetReferences(refCount, refImg);

    zetaFilt.SetPatchRadius(patchRadius);

    zetaFilt.SetNeighbourhoodRadius(nbhdRadius);

    zetaFilt.Initialise();


}




