#ifndef CAMERA_H
#define CAMERA_H

#include "ray.h"

class Camera 
{
public:

	Vector3d mOrigin;
	Vector3d mImagePlaneTopLeftCorner;
	Vector3d mImagePlaneWidth;;
	Vector3d mImagePlaneHeight;

	Camera( Vector3d pImagePlaneTopLeftCorner, Vector3d pImagePlaneWidth, Vector3d pImagePlaneHeight)
	{
		mImagePlaneHeight = pImagePlaneHeight;
		mImagePlaneWidth = pImagePlaneWidth;
		mImagePlaneTopLeftCorner = pImagePlaneTopLeftCorner;
	}

	Ray getRayP1(float pPixPosX, float pPixPosY,Vector3d pDirection)
	{
		mOrigin = mImagePlaneTopLeftCorner + (pPixPosX * mImagePlaneWidth) + (pPixPosY * mImagePlaneHeight);
		Ray lRay(mOrigin, pDirection);
		return lRay;
	}
	
	Ray getRayP3(float pPixPosX, float pPixPosY, Vector3d pOrigin)
	{
		Vector3d lDirection = mImagePlaneTopLeftCorner + (pPixPosX * mImagePlaneWidth) + (pPixPosY * mImagePlaneHeight);
		Ray lRay(pOrigin, lDirection);
		return lRay;
	}
	
};
#endif // !CAMERA_H
