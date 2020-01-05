#ifndef SPHERE_H
#define SPHERE_H

#include "shapes.h"
#include "objtype.h"

class Sphere : public Shapes
{

public:
	Vector3d mCenter;
	double mRadius;
	Vector3d mColor;
	ObjType mObjType;

	Sphere(Vector3d pCenter, double pRadius, Vector3d pColor)
	{
		mCenter = pCenter;
		mRadius = pRadius;
		mColor = pColor;
	}

	Sphere(Vector3d pCenter, float pRadius, Vector3d pColor, ObjType pObjType)
	{
		mCenter = pCenter;
		mRadius = pRadius;
		mObjType = pObjType;
		mColor = pColor;
	}

	bool intersect(Ray &pRay, double pTMin, double pTMax, HitCount &pHitCount) const
	{
		Vector3d lOriginToCenter = pRay.getOrigin() - mCenter;

		double lA = pRay.getDirection().dot(pRay.getDirection());

		double lB =  2 * lOriginToCenter.dot(pRay.getDirection());

		double lC = ((lOriginToCenter.dot(lOriginToCenter)) - ( mRadius * mRadius ));;

		double lDiscriminant = ( (lB * lB) - ( 4 * (lA * lC) ) ) ;


		if (lDiscriminant > 0)
		{
			double lT1 = (-lB - sqrt((lB * lB) - (4 * (lA * lC)))) / (2 * lA);

			if ( lT1 < pTMax && lT1 > pTMin )
			{
				pHitCount.setT( lT1 );
				pHitCount.setIntersectionPoint( pRay.AtPoint( lT1 ));
				pHitCount.setNormal( ( ( pHitCount.getIntersectionPoint() - mCenter ) / mRadius ) );
				pHitCount.mObjType = mObjType;
				pHitCount.setObjectColor(mColor);
				return true;
			}
			double lT2 = (-lB + sqrt((lB * lB) - (4 * (lA * lC)))) / (2 * lA);
			if ( lT2 < pTMax && lT2 > pTMin)
			{
				pHitCount.setT( lT2 );
				pHitCount.setIntersectionPoint( pRay.AtPoint( lT2 ));
				pHitCount.setNormal( ( ( pHitCount.getIntersectionPoint() - mCenter ) / mRadius ) );
				pHitCount.mObjType = mObjType;
				pHitCount.setObjectColor( mColor );
				return true;
			}
		}
		return false;
	}


};

#endif // !SPHERE_H
