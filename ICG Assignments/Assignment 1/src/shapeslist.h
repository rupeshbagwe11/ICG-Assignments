#ifndef SHAPESLIST_H
#define SHAPESLIST_H

#include "shapes.h"

class ShapesList : public Shapes
{
public:
	Shapes** mShapesList;
	int miListSize;

	ShapesList(Shapes** pShapesList, int pSize )
	{
		mShapesList = pShapesList;
		miListSize = pSize;
	}

	bool intersect(Ray &pRay, double pTMin, double pTMax, HitCount &pHitCount) const
	{
		HitCount lTemp;
		bool lHit = false;
		double closest = pTMax;

		for (int i = 0; i < miListSize; i++)
		{
			if (mShapesList[i]->intersect(pRay, pTMin, closest, lTemp))
			{
				lHit = true;
				closest = lTemp.getT();
				pHitCount = lTemp;
			}
		}
		return lHit;

	}


};

#endif // !WORLD_H