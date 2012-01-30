/*
 * abcdRandomRotationSource.cc
 *
 *  Created on: Dec 2, 2011
 *      Author: paul
 */


#include <sys/types.h>


#ifndef WIN32
#include <sys/time.h>
#endif

#include <abcdRandomRotationSource.h>

abcdRandomRotationSource::abcdRandomRotationSource()
{
	this->_homogeneous = false;

#ifndef WIN32
  timeval tv;
  gettimeofday(&tv, NULL);
  this->_randInit = tv.tv_usec;
#else
  cerr << "abcdRandomRotationSource::abcdRandomRotationSource: Not implemented yet for Windows" << endl;
  exit(1);
#endif

}


abcdRandomRotationSource::~abcdRandomRotationSource()
{

}

irtkMatrix abcdRandomRotationSource::Get()
{
	irtkMatrix R(3,3);
	irtkVector u2(3), u1(3), u0(3), w(3);

	int i,j;
	double phi, costheta, sintheta, theta_w;

	// Spherical coordinates for u2 random point on unit sphere's surface.
	// Angle from x-axis within x-y plane
	phi = 2.0f * M_PI * ran2(& this->_randInit);
	// Angle from north pole (z=1)
	costheta = 2 * ran2(& this->_randInit) - 1.0;
	sintheta = sqrt(1 - costheta*costheta);

	// Cartesian coordinates of u2
	u2(0) = cos(phi) * sintheta;
	u2(1) = sin(phi) * sintheta;
	u2(2) = costheta;

	// Obtain a starting vector w that is orthogonal to u2 and
	// in the x-y plane
	if (abs(u2(0)) < FLT_MIN){
		w(0) = 1;
		w(1) = 0;
		w(2) = 0;
	} else {
		w(0) = u2(1);
		w(1) = -1.0*u2(0);
		w(2) = 0;
		w.Normalize();
	}

	// Randomly chosen angle of a rotation around u2 to apply to w
	theta_w = 2.0f * M_PI * ran2(& this->_randInit);

	irtkMatrix outerProduct(3,3);
	irtkMatrix skewSym(3,3);
	irtkMatrix I(3,3);

	I.Ident();

	// Reset skew symmetric matrix.
	skewSym = I - I;
	// And assign entries so that the skew symmetric matrix has the the same effect as
	// the cross product u2 x v when premultiplied by any vector v, i.e
	// skewSym * v == u x v.
	skewSym(0, 2) = u2(1);
	skewSym(1, 0) = u2(2);
	skewSym(2, 1) = u2(0);

	skewSym(2, 0) = -1.0 * u2(1);
	skewSym(0, 1) = -1.0 * u2(2);
	skewSym(1, 2) = -1.0 * u2(0);

	// Outer product u2 * u2'
	for (j = 0; j < 3; j++){
		for (i = 0; i < 3; i++){
			outerProduct(i,j) = u2(i) * u2(j);
		}
	}

	// Matrix form for a rotation of theta_w around vector u2
	R = outerProduct +  ((I - outerProduct) * cos(theta_w)) + skewSym * sin(theta_w);

	// Rotate w with this to get u1
	u1 = R * w;


	// Now we can get u0
	u0 = u1.CrossProduct(u2);

	for (i = 0; i < 3; i++){
		R(i,0) = u0(i);
		R(i,1) = u1(i);
		R(i,2) = u2(i);
	}

	irtkMatrix temp;


//	cerr << "===================" << endl;
//	u2.Print();
//	R.Print();
//
//	temp = R;
//	temp.Transpose();
//	temp = temp * R;
//	temp.Print();
//
//	temp = R;
//	temp.Transpose();
//	temp = R * temp;
//	temp.Print();
//	cerr << "===================" << endl;

	if (_homogeneous){
		irtkMatrix Rh(4,4);
		Rh(R, 0,0);
		Rh(3,3) = 1;
		return Rh;
	}

	return R;

}
