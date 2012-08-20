/*
 * abcdPointDistanceMap.cc
 *
 *  Created on: Aug 3, 2012
 *      Author: paul
 */


#include <irtkImage.h>

#include <vector>
#include <abcdPointDistanceMap.h>

template <class VoxelType> abcdPointDistanceMap<VoxelType>::abcdPointDistanceMap()
{

  this->_numberOfSeedPoints = 0;

}


template <class VoxelType> abcdPointDistanceMap<VoxelType>::~abcdPointDistanceMap(void)
{

}

template <class VoxelType> bool abcdPointDistanceMap<VoxelType>::RequiresBuffering(void)
{
  return true;
}


template <class VoxelType> const char *abcdPointDistanceMap<VoxelType>::NameOfClass()
{
  return "abcdPointDistanceMap";
}

template <class VoxelType> void abcdPointDistanceMap<VoxelType>::Initialize()
{
	// Do the initial set up
	this->irtkImageToImage<VoxelType>::Initialize();

  irtkConnectivityType conn = CONNECTIVITY_06;
  this->_faceOffsets.Initialize(this->_input, conn);

  if (this->_input->GetT() > 1){
	  cerr << "Warning: abcdPointDistanceMap<VoxelType>::Initialize : Image is 4D, ignoring all volumes after first." << endl;
	}
}

template <class VoxelType> void abcdPointDistanceMap<VoxelType>::Finalize()
{
  // Base class:
  this->irtkImageToImage<VoxelType>::Finalize();

}

/// Add a point in world coordinates
template <class VoxelType> void abcdPointDistanceMap<VoxelType>::addSeedPointW(double x, double y, double z)
{
	if (this->_numberOfSeedPoints >= MAX_POINTS)
	{
		cerr << "abcdPointDistanceMap<VoxelType>::addSeedPointW : Point limit exceeded." << endl;
		exit(1);
	}

	this->_input->WorldToImage(x, y, z);
	this->addSeedPointI( (int) x, (int) y, (int) z);

}

/// Add a point in image coordinates
template <class VoxelType> void abcdPointDistanceMap<VoxelType>::addSeedPointI(int i, int j, int k)
{
	if (this->_numberOfSeedPoints >= MAX_POINTS)
	{
		cerr << "abcdPointDistanceMap<VoxelType>::addSeedPointI : Point limit exceeded." << endl;
		exit(1);
	}

	_seedpoints_i[_numberOfSeedPoints] = (int) i;
	_seedpoints_j[_numberOfSeedPoints] = (int) j;
	_seedpoints_k[_numberOfSeedPoints] = (int) k;

	++_numberOfSeedPoints;


}

template <class VoxelType> void abcdPointDistanceMap<VoxelType>::Run()
{
  int i, j, k, n, dimx, dimy, dimz, currVal, currBoundarySize, faceNbrInd;
  VoxelType *ptr2output, *ptr2input;
  bool modified;


  this->Initialize();

  ptr2output = this->_output->GetPointerToVoxels();

  // Clear output.
  for (n = 0; n < this->_output->GetNumberOfVoxels(); ++n){
    *ptr2output = 0;
    ++ptr2output;
  }

  if (this->_numberOfSeedPoints == 0){
	  this->Finalize();
		return;
	}

  dimx = this->_input->GetX();
  dimy = this->_input->GetY();
  dimz = this->_input->GetZ();

  // Clear any voxels at boundary of field of view to prevent any boundary errors.
  for (i = 0; i < dimx; ++i){
    for (j = 0; j < dimy; ++j){
      this->_input->Put(i, j, 0       , 0);
      this->_input->Put(i, j, dimz - 1, 0);
    }
  }
  for (j = 0; j < dimy; ++j){
    for (k = 0; k < dimz; ++k){
      this->_input->Put(0       , j, k, 0);
      this->_input->Put(dimx - 1, j, k, 0);
    }
  }
  for (i = 0; i < dimx; ++i){
    for (k = 0; k < dimz; ++k){
      this->_input->Put(i, 0       , k, 0);
      this->_input->Put(i, dimy - 1, k, 0);
    }
  }


  vector<int> currVoxels(this->_numberOfSeedPoints, 0);
  vector<int> nextVoxels;
  vector<int>::iterator last;

  currVal = 1;

  for (n = 0; n < this->_numberOfSeedPoints; ++n){
    currVoxels[n] = this->_input->VoxelToIndex(_seedpoints_i[n], _seedpoints_j[n], _seedpoints_k[n]);
    this->_output->Put(_seedpoints_i[n], _seedpoints_j[n], _seedpoints_k[n], currVal);

  }

  ptr2output = this->_output->GetPointerToVoxels();
  ptr2input  = this->_input->GetPointerToVoxels();



  do {
    // Loop over all voxels in the current boundary.  Look at each of its
    // neighbours. If it is in the mask, mark it with the next distance
    // value.  The largest number of voxels for the next round is achieved
    // when all voxels in the current round have completely separate face
    // neighbours.

    // Likely not to need all this storage.  There should be duplicates,
    // i.e. face neighbours shared by more than one of the current boundary
    // voxels.

    modified = false;

    currBoundarySize = (int) currVoxels.size();

    nextVoxels.assign(6 * currBoundarySize , 0);

    for (n = 0; n < currBoundarySize; ++n){
      for (i = 0; i < 6; ++i){
        faceNbrInd = currVoxels[n] + this->_faceOffsets(i);

        if (( *(ptr2input  + faceNbrInd) > 0) &&
            ( *(ptr2output + faceNbrInd) < 1) ){
          modified = true;
          nextVoxels[6*n + i] = faceNbrInd;
        }

      } // face neighbour loop.
    } // boundary voxel loop.

    std::sort(nextVoxels.begin(), nextVoxels.end());
    // Puts the unique values at the front of the vector, potential garbage in remaining positions.
    last = std::unique(nextVoxels.begin(), nextVoxels.end());
    nextVoxels.erase(last, nextVoxels.end());

    if (modified){

      currVoxels = nextVoxels;
      ++currVal;

      currBoundarySize = (int) currVoxels.size();

      for (n = 0; n < currBoundarySize; ++n){
        *(ptr2output + currVoxels[n]) = currVal;
      }

    }

  } while (modified);

  this->_output->Write("bla.nii.gz");

}


template class abcdPointDistanceMap<irtkGreyPixel>;






