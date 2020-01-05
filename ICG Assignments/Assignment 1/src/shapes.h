#ifndef SHAPES_H
#define SHAPES_H

#include "ray.h"
//#include "objecttype.h"
#include "hitcount.h"

class Shapes
{
public:
	virtual bool intersect(Ray &pRay, double pTMin, double pTMax, HitCount &pHitCount) const = 0;
};
#endif // !SHAPES_H
