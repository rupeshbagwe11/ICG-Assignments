#ifndef RAY_H
#define RAY_H

#include "utils.h"

using namespace Eigen;

class Ray
{

private:
	Vector3d mOrigin;
	Vector3d mDirection;

public:
	Ray() {}
	Ray(Vector3d pOrigin, Vector3d pDirection)
	{
		mOrigin = pOrigin;
		mDirection = pDirection;
	}

	Vector3d getOrigin()
	{
		return mOrigin;
	}

	Vector3d getDirection()
	{
		return mDirection;
	}

	Vector3d AtPoint(float pT)
	{
		return ( mOrigin + (pT * mDirection) );
	}

};

#endif // !RAY_H
