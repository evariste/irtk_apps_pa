/*
 * abcdRandomRotationSource.h
 *
 *  Created on: Dec 2, 2011
 *      Author: paul
 */

// Julie C Mitchell, Sampling Rotation Groups by Successive Orthogonal
// Image, SIAM J Sci comput. 30(1), 2008, pp 525-547

// Seek R = [u0 u1 u2] where the columns are the images of the unit axis
// vectors under the rotation.

// u2 uniformly sampled from a sphere (see
// http://mathworld.wolfram.com/SpherePointPicking.html)


#ifndef ABCDRANDOMROTATIONSOURCE_H_
#define ABCDRANDOMROTATIONSOURCE_H_

#include <irtkImage.h>
#include <irtkVector.h>
#include <irtkMatrix.h>

#include <gsl/gsl_rng.h>
//#include <nr.h>


class abcdRandomRotationSource : public irtkObject
{

	// Seed for random number generator.
	long int _randInit;
	
	gsl_rng * _r; 
	const gsl_rng_type * _T;

	// Require a 4 by 4 homogeneous matrix if true
	bool _homogeneous;


public:

	abcdRandomRotationSource();

	virtual ~abcdRandomRotationSource();

	// Get a random rotation matrix;
	irtkMatrix Get();

	// Set the homogenous flag
	void SetHomogeneous(bool);

};


inline void abcdRandomRotationSource::SetHomogeneous(bool val){

	this->_homogeneous = val;

}


#endif /* ABCDRANDOMROTATIONSOURCE_H_ */
