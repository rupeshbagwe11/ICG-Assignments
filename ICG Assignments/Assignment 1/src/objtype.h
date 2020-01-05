#ifndef OBJTYPE_H
#define OBJTYPE_H

class ObjType
{
public:
	float mDiffuse;
	float mAmbient;
	float mSpecular;
	ObjType()
	{
		mDiffuse = 1.0;
		mAmbient = 1.0;
		mSpecular = 1.0;
	}
	ObjType(float pDiffuse, float pAmbient, float pSpecular)
	{
		mDiffuse = pDiffuse;
		mAmbient = pAmbient;
		mSpecular = pSpecular;
	}
};

#endif // !1



