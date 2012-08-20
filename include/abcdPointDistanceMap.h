/*
 * abcdPointDistanceMap.h
 *
 *  Created on: Aug 3, 2012
 *      Author: paul
 */

#ifndef ABCDPOINTDISTANCEMAP_H_
#define ABCDPOINTDISTANCEMAP_H_

#include <irtkImageToImage.h>

#define MAX_POINTS 100


template <class VoxelType> class abcdPointDistanceMap : public irtkImageToImage <VoxelType>{

protected:

	/// Image grid locations for starting points
	int _seedpoints_i[MAX_POINTS];
	int _seedpoints_j[MAX_POINTS];
	int _seedpoints_k[MAX_POINTS];

	int _numberOfSeedPoints;

	irtkNeighbourhoodOffsets _faceOffsets;

	/// Returns the name of the class
	virtual const char *NameOfClass();

  /// Initialize the filter
  virtual void Initialize();

  /// Finalize the filter
  virtual void Finalize();

	/// Returns whether the filter requires buffering
	virtual bool RequiresBuffering();

public:
	/// Constructor
	abcdPointDistanceMap();

	/// Destructor
	virtual ~abcdPointDistanceMap();

	/// Run filter
	virtual void Run();


	/// Add a point in world coordinates
	void addSeedPointW(double, double, double);

	/// Add a point in image coordinates
	void addSeedPointI(int, int, int);


};

#endif /* ABCDPOINTDISTANCEMAP_H_ */
