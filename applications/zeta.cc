

#include "irtkImage.h"

#include "zeta.h"

#ifdef HAS_MPI
#include <mpi.h>
#endif


char *input_name = NULL, *output_name = NULL;
char *mask_name = NULL;
char **ref_name;




void usage()
{
    cout << " " << endl;
    cout << " Obtain a zeta map for a target image based on a given set of reference images." << endl;
    cout << "Usage" << endl;
    cout << " " << endl;
    cout << " zeta [Target] [Output] [-refs Ref1 Ref2 ... RefN] <options>" << endl;
    cout << " " << endl;
    cout << " [Target]        : Target or query image." << endl;
    cout << " [Output]        : Where to save the zeta map." << endl;
    cout << " Ref1, ..., RefN : N reference images (if 4D volumes then channels must " << endl;
    cout << "                   match channels in target)." << endl;
    cout << " Options:" << endl;
    cout << " " << endl;
    cout << " -patchRadius Rp : Radius of patch (integer). p=1 gives a 3x3x3 patch size, " << endl;
    cout << "                   p=2 gives 5x5x5 etc." << endl;
    cout << " -nbhdRadius Rn  : Radius of neighbourhood (integer) over which to search " << endl;
    cout << "                   patches in the reference set." << endl;
    cout << " -k K            : Number of nearest neighbours to use for zeta calculation " << endl;
    cout << "                   (integer)  K <= N (see above)." << endl;
    cout << " -mask image     : Mask image to define region of interest x,y,z dimensions " << endl;
    cout << "                   match those in target." << endl;
    cout << " -mahalanobis    : Use Mahalanobis distance between target and reference " << endl;
    cout << "                   patches (this is the default so flag does not really " << endl;
    cout << "                   need to be set)." << endl;
    cout << " -euclidean      : Use Euclidean distance between target and reference " << endl;
    cout << "                   patches. Overrides default Mahalanobis distance." << endl;

    exit(0);
}


int main(int argc, char **argv)
{

    int ok;

    int patchRadius = 0; // Diam = 1 + 2*R
    int nbhdRadius = 4;
    int refCount = 0;

    int k = 3;

    bool MahalanobisFlag = true;

#ifdef HAS_MPI
    /* Initialize */
    MPI_Init(&argc,&argv);
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

    if (argc < 2){
#ifdef HAS_MPI
      MPI_Finalize();
#endif
        usage();
    }


    input_name  = argv[1];
    argc--;
    argv++;
    output_name = argv[1];
    argc--;
    argv++;

#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cout << "Target image: " << input_name << endl;
    cout << "Saving output to: " << output_name << endl;
    }
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
#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cout << "Setting patch size to " << patchRadius << endl;
        }
        ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-nbhdRadius") == 0)) {
        argc--;
        argv++;
        nbhdRadius = atoi(argv[1]);
        argc--;
        argv++;

#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cout << "Setting neighbourhood size to " << nbhdRadius << endl;
        }
        ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-k") == 0)) {
        argc--;
        argv++;
        k = atoi(argv[1]);
        argc--;
        argv++;
#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cout << "Setting number of neighbours (k) to " << k << endl;
        }
        ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
        argc--;
        argv++;
        mask_name = argv[1];
        argc--;
        argv++;
#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cout << "Using mask file: " << mask_name << endl;
        }
        ok = true;
      }

      if ((ok == false) && (strcmp(argv[1], "-euclidean") == 0)) {
        argc--;
        argv++;
        MahalanobisFlag = false;
#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cout << "Using Euclidean distance." << endl;
        }
        ok = true;
      }

      if ((ok == false) && (strcmp(argv[1], "-mahalanobis") == 0)) {
        argc--;
        argv++;
        MahalanobisFlag = true;
#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cout << "Using Mahalanobis distance." << endl;
        }
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
#ifdef HAS_MPI
        if (myid == 0)
#endif
        {
        cerr << "Unknown option: " << argv[1] << endl;
        }
        usage();
      }
    }


    if (refCount == 0){
#ifdef HAS_MPI
      if (myid == 0)
#endif
      {
      cout << "No reference images. Nothing to do." << endl;
      cout << "Exiting" << endl;
      }
      exit(0);
    }

    irtkRealImage **refImg = new irtkRealImage*[refCount];
#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    cout << "Using the reference images: " << endl;
    for (int i = 0; i < refCount; i++){
      cout << "     " << ref_name[i] << endl;
    }
    }

    for (int i = 0; i < refCount; i++){
      refImg[i] = new irtkRealImage(ref_name[i]);
    }

    Zeta zetaFilt;

    zetaFilt.SetTarget(&image);
    zetaFilt.SetReferences(refCount, refImg);
    zetaFilt.SetPatchRadius(patchRadius);
    zetaFilt.SetNeighbourhoodRadius(nbhdRadius);
    zetaFilt.SetK(k);
    zetaFilt.UseMahalanobis(MahalanobisFlag);

    if (mask_name != NULL){
      irtkGreyImage *mask = new irtkGreyImage(mask_name);
      zetaFilt.SetMask(mask);
    }

    zetaFilt.Initialise();

#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    zetaFilt.Print();
    }

#ifdef HAS_MPI
    zetaFilt.RunParallel();
#else
    zetaFilt.Run();
#endif

#ifdef HAS_MPI
    if (myid == 0)
#endif
    {
    irtkRealImage *out = zetaFilt.GetOutput();
    out->Write(output_name);
    }

#ifdef HAS_MPI
    MPI_Finalize();
#endif


}




