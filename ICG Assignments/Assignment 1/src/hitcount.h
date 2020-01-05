#ifndef HITCOUNT_H
#define HITCOUNT_H

#include "utils.h"
#include "objtype.h"

using namespace Eigen;



class HitCount
{
public:
	float mT; //Intersection Point
	Vector3d mP; //Vector at Point T
	Vector3d mNormal; //Vector Normal to Vector P
	ObjType mObjType;
	Vector3d mObjectColor;

	HitCount()
	{
	}

	float getT()
	{
		return mT;
	}
	void setT(float pT)
	{
		mT = pT ;
	}
	Vector3d getIntersectionPoint()
	{
		return mP;
	}
	void setIntersectionPoint(Vector3d pP)
	{
		mP = pP;
	}
	Vector3d getNormal()
	{
		return mNormal;
	}
	void setNormal(Vector3d pNormal)
	{
		mNormal = pNormal;
	}
	Vector3d getObjectColor()
	{
		return mObjectColor;
	}
	void setObjectColor(Vector3d pObjectColor)
	{
		mObjectColor = pObjectColor;
	}
	ObjType getObjectType()
	{
		return mObjType;
	}
	void setObjectType(ObjType pObjectType)
	{
		mObjType = pObjectType;
	}

};
#endif // !HITCOUNT_H
