#include <irtkImage.h>
//#include <nr.h>
#include <gsl/gsl_rng.h>

#include <sys/types.h>


//#ifdef WIN32
//#include <time.h>
//#else
//#include <sys/time.h>
//#endif


#include <map>



char *output_name = NULL, **input_names = NULL;
char **weight_names = NULL;
char *mask_name = NULL;

//long ran2Seed;
//long ran2initialSeed;
gsl_rng * ranGen; 
const gsl_rng_type * ranGenType;

// first short is the label, second float is the combined weight from images votes
// for that label.
typedef map<short, float> countMap;

map<short, float>::iterator iter;

short getMostPopular(countMap cMap){
  double maxWeight = 0, mostPopLabel = -1;

  if (cMap.size() == 0){
    // No votes to count, treat as background.
    return 0;
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second > maxWeight){
      maxWeight     = iter->second;
      mostPopLabel = iter->first;
    }
  }

  return mostPopLabel;
}

bool isEquivocal(countMap cMap){
  double maxWeight = 0;
  short numberWithMax = 0;

  if (cMap.size() == 0){
    // No votes to count, treat as background.
    return false;
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second > maxWeight){
      maxWeight     = iter->second;
    }
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (fabs(iter->second - maxWeight) < 1E-09){
      ++numberWithMax;
    }
  }

  if (numberWithMax > 1){
    return true;
  } else {
    return false;
  }

}

short decideOnTie(countMap cMap){
  double maxWeight = 0;
  short numberWithMax = 0;
  int index, count;
  short temp;
  double val;

  if (cMap.size() == 0){
    // No votes to count, treat as background.
    return false;
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second > maxWeight){
      maxWeight     = iter->second;
    }
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (fabs(iter->second -  maxWeight) < 1E-09){
      ++numberWithMax;
    }
  }

  short *tiedLabels = new short[numberWithMax];

  count = 0;
  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (fabs(iter->second -  maxWeight) < 1E-09){
      tiedLabels[count] = iter->first;
      ++count;
    }
  }

  //val =  ran2(&ran2Seed);
  val =  gsl_rng_uniform(ranGen);
  index = (int) round(val * (count - 1));

  temp = tiedLabels[index];

  delete tiedLabels;

  return temp;
}


void usage()
{
  cerr << "Usage: combineLabels2 [output] [N] [input1] .. [inputN] [weight1] ... [weightN] <-options>" << endl;
  cerr << "The images [input1] .. [inputN] are assumed to contain labellings that" << endl;
  cerr << "are in register.  For each voxel, the modal label from the images is used " << endl;
  cerr << "to generate a label in the output image." << endl;
  cerr << "-pad [value]  : Padding value, default = -1." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int numberOfClassifiers, i, v, voxels, contendedVoxelCount, contendedVoxelIndex, equivocalCount;
  bool ok;

  irtkGreyPixel *pIn, *pIn_0, *pOut;
  //  irtkGreyPixel minLabel = MAX_GREY, maxLabel = MIN_GREY;
  irtkGreyImage input, output;
  irtkGreyImage *input_0;
  irtkByteImage mask;
  irtkBytePixel *pMask;
  bool writeMask = false;
  int pad = -1;
  irtkRealImage weightImg;
  irtkRealPixel *ptr2weight;

  // Check command line
  if (argc < 4){
    usage();
  }

  output_name = argv[1];
  argc--;
  argv++;

  numberOfClassifiers = atoi(argv[1]);
  argc--;
  argv++;

  input_names = new char *[numberOfClassifiers];
  for (i = 0; i < numberOfClassifiers; i++){
    input_names[i] = argv[1];
    argc--;
    argv++;
  }

  weight_names = new char *[numberOfClassifiers];
  for (i = 0; i < numberOfClassifiers; i++){
      weight_names[i] = argv[1];
      argc--;
      argv++;
  }

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;      argv++;
      pad = atoi(argv[1]);
      argc--;      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Mask has a value of 1 for all uncontended voxels.
  mask.Read(input_names[0]);
  voxels = mask.GetNumberOfVoxels();

  // Reset mask so that all voxels are marked as unanimous.
  pMask = mask.GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    *pMask = 1;
    ++pMask;
  }


  input_0 = new irtkGreyImage(input_names[0]);

  for (i = 1; i < numberOfClassifiers; ++i){

    input.Read(input_names[i]);

    pIn   = input.GetPointerToVoxels();
    pIn_0 = input_0->GetPointerToVoxels();
    pMask = mask.GetPointerToVoxels();

      for (v = 0; v < voxels; ++v){

        if (*pIn != *pIn_0 && (*pIn > pad || *pIn_0 > pad) ){
          *pMask = 0;
        }
        ++pIn;
        ++pIn_0;
        ++pMask;
      }
  }

  // free some memory.
  delete input_0;

  contendedVoxelCount = 0;
  pMask = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pMask == 0){
      ++contendedVoxelCount;
    }
    ++pMask;
  }

  cout << "Number of contended Voxels = " << contendedVoxelCount << " = ";
  cout << 100.0 * contendedVoxelCount / ((double) voxels) << "%" << endl;

  countMap *counts = new countMap[contendedVoxelCount];

  for (i = 0; i < numberOfClassifiers; ++i){

    input.Read(input_names[i]);
    weightImg.Read(weight_names[i]);

	pIn = input.GetPointerToVoxels();
    pMask = mask.GetPointerToVoxels();
    ptr2weight = weightImg.GetPointerToVoxels();

    contendedVoxelIndex = 0;
    for (v = 0; v < voxels; ++v){

      //Is the current voxel contentded?
      if (*pMask == 0){
        if (*pIn > pad){
          counts[contendedVoxelIndex][*pIn] += *ptr2weight;
        }
        ++contendedVoxelIndex;
      }
      ++pIn;
      ++pMask;
      ++ptr2weight;
    }
  }

  // Initialise output image using first input image, all inputs assumed
  // same size etc..
  output.Read(input_names[0]);
  pOut  = output.GetPointerToVoxels();
  pMask = mask.GetPointerToVoxels();

  contendedVoxelIndex = 0;
  for (v = 0; v < voxels; ++v){

    //Is the current voxel contentded?
    if (*pMask == 0){
      *pOut = getMostPopular(counts[contendedVoxelIndex]);
      ++contendedVoxelIndex;
    }

    ++pOut;
    ++pMask;
  }

  // Get ready for random stuff.
  //time_t tv;
  //tv = time(NULL);

  //ran2Seed = tv; // .tv_usec;
  //ran2initialSeed = -1 * ran2Seed;
  //(void) ran2(&ran2initialSeed);

  gsl_rng_env_setup();
  ranGenType = gsl_rng_default;
  ranGen = gsl_rng_alloc (ranGenType);


  contendedVoxelIndex = 0;
  equivocalCount = 0;

  pOut  = output.GetPointerToVoxels();
  pMask = mask.GetPointerToVoxels();
  for (v = 0; v < voxels; ++v){

    //Is the current voxel contentded?
    if (*pMask == 0){
      if (isEquivocal(counts[contendedVoxelIndex])){
        ++equivocalCount;
        *pOut = decideOnTie(counts[contendedVoxelIndex]);
      }
      ++contendedVoxelIndex;
    }
    ++pOut;
    ++pMask;
  }

  cout << "Number of equivocal voxels = " << equivocalCount << " = ";
  cout << 100.0 * equivocalCount / ((double) voxels) << "%" << endl;

  output.Write(output_name);

  delete [] counts;
}
