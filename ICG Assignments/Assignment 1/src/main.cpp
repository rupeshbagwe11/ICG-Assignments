// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <fstream>


// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include "ray.h"
#include "camera.h"
#include "sphere.h"
#include "shapeslist.h"
#include "objtype.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;



Vector3d getColorP1(Ray pRay, Shapes* pWorld, Vector3d pLightSource)
{
	Vector3d WhiteColor(1.0, 1.0, 1.0);
	Vector3d GreenColor(0.2, 1.0, 0.2);
	Vector3d RedColor(1.0, 0.0, 0.0);
	Vector3d BlueColor(0.0, 0.0, 1.0);

	double lMaxFloatNumber = getMaxFloatNumber();

	HitCount lHitCount;
	if ( pWorld -> intersect ( pRay, 0.001, lMaxFloatNumber, lHitCount ) )
	{
		Vector3d lRayIntersection = lHitCount.getIntersectionPoint();
		Vector3d lRayNormal = lRayIntersection.normalized();
		// Simple diffuse model
		double lcolor = (pLightSource - lRayIntersection).normalized().transpose() * lRayNormal;
		// Clamp to zero
		lcolor = max(lcolor, 0.);
		if (lcolor > 1)
		{
			lcolor = 1.0;
		}

		return (lcolor * lHitCount.getObjectColor());
	}
	else {
		return BlueColor;
	}

}



void part1()
{
	std::cout << "Ray Tracing Spheres" << std::endl;

	const std::string filename("part1.png");

	int liWidth = 800;
	int liHeight = 800;
	//int liAvg = 50;

	MatrixXd lColorRedValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorGreenValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorBlueValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lAlphaValue = MatrixXd::Zero(liWidth, liHeight); // Store the alpha mask

	Vector3d lImagePlaneWidth(4.0, 0.0, 0.0);
	Vector3d lImagePlaneHeight(0.0, -4.0, 0.0);
	Vector3d lTopLeftCorner(-2.0, 2.0, 0);
	Vector3d lDirection(0, 0, -1.0);

	Camera lCamera( lTopLeftCorner, lImagePlaneWidth, lImagePlaneHeight);

	Vector3d RedColor(1.0, 0.0, 0.0);
	Vector3d BlueColor(0.0, 0.0, 1.0);
	Vector3d GreenColor(0.0, 1.0, 0.0);
	Vector3d YellowColor(1.0, 1.0, 0.4);
	
	Shapes* lObjlist[3];

	Vector3d lSphere1Center(0.5, 1.0, 0.0);
	Vector3d lSphere2Center(-1.2, 0.9, 0.0);
	Vector3d lSphere3Center(-0.4, 0.0, 0.0);

	lObjlist[0] = new Sphere(lSphere1Center, 0.5, RedColor);
	lObjlist[1] = new Sphere(lSphere2Center, 0.3, GreenColor);
	lObjlist[2] = new Sphere(lSphere3Center, 0.7, YellowColor);

	Shapes* lWorld = new ShapesList(lObjlist, 3);
	Vector3d lLightSource(-1.0, 2.0, -1.0);

	for (int i = 0; i < liWidth; i++)
	{
		for (int j = 0; j < liHeight; j++)
		{

		//	Vector3d lColor(0, 0, 0);

			float lfPixPosX = float(i) / float(liWidth);
			float lfPixPosY = float(j) / float(liHeight);

			Ray lRay = lCamera.getRayP1(lfPixPosX, lfPixPosY, lDirection);

		//	for (int k = 0; k < liAvg; k++)
		//	{
		//		lColor = lColor + getColorP1(lRay, lWorld);
		//	}

		//	Vector3d lFColor = lColor / liAvg;
			Vector3d lFColor = getColorP1(lRay, lWorld, lLightSource);

			//Vector3d lFColor = lColor * RedColor;

			lColorRedValue(i, j) = lFColor[0];
			lColorGreenValue(i, j) = lFColor[1];
			lColorBlueValue(i, j) = lFColor[2];

			// Disable the alpha mask for this pixel
			lAlphaValue(i, j) = 1;
		}

	}


	// Write it in a png image. Note that the alpha channel is reversed to make the white (color = 1) pixels transparent (alhpa = 0)
	write_matrix_to_png(lColorRedValue, lColorGreenValue, lColorBlueValue, lAlphaValue, filename);

}


double mdAmbient = 0.5;
int miPhongsExponent = 16;


Vector3d getColorP2(Ray pRay, Shapes* pWorld, Vector3d pLightSources[], int pLSsize, Vector3d pEyeSource)
{
	Vector3d WhiteColor(1.0, 1.0, 1.0);
	Vector3d BlueColor(0.0, 0.0, 1.0);

	double lMaxFloatNumber = getMaxFloatNumber();

	HitCount lHitCount;
	if (pWorld->intersect(pRay, 0.001, lMaxFloatNumber, lHitCount))
	{
		double lcolor = 0.0;
		for (int m = 0; m < pLSsize; m++)
		{

			Vector3d lRayIntersection = lHitCount.getIntersectionPoint();
			Vector3d lRayNormal = lRayIntersection.normalized();

			ObjType lObjType = lHitCount.getObjectType();

			//------diffuse-------
			// Simple diffuse model
			double lColorByDiffuse = (pLightSources[m] - lRayIntersection).normalized().transpose() * lRayNormal;

			// Clamp to zero
			double lFColorByDiffuse = lObjType.mDiffuse * clampNumberBetween0and1(lColorByDiffuse);

			//------specular-------
			Vector3d lV = pEyeSource - lRayIntersection;
			Vector3d lL = pLightSources[m] - lRayIntersection;
			Vector3d lH = (lV + lL).normalized();
			double lColorBySpecular = lH.normalized().transpose() * lRayNormal;
			double lFColorBySpecular = lObjType.mSpecular * (pow(clampNumberBetween0and1(lColorBySpecular), miPhongsExponent));

			//------ambient-------
			double lFColorByAmbient = lObjType.mAmbient * mdAmbient;

			//-----------------
			double lColorByAllThree = lFColorByAmbient + lFColorByDiffuse + lFColorBySpecular;
			double lFColorByAllThree = clampNumberBetween0and1(lColorByAllThree);

			lcolor = lcolor + lFColorByAllThree;
		}

		Vector3d lFColor = clampNumberBetween0and1(lcolor) * lHitCount.getObjectColor();


		return lFColor;

	}
	else {
		return BlueColor;
	}

}




void part2()
{
	std::cout << "Shading" << std::endl;

	const std::string filename("part2.png");

	int liWidth = 800;
	int liHeight = 800;
	//int liAvg = 50;

	MatrixXd lColorRedValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorGreenValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorBlueValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lAlphaValue = MatrixXd::Zero(liWidth, liHeight); // Store the alpha mask

	Vector3d lImagePlaneWidth(4.0, 0.0, 0.0);
	Vector3d lImagePlaneHeight(0.0, -4.0, 0.0);
	Vector3d lTopLeftCorner(-2.0, 2.0, 0);
	Vector3d lDirection(0, 0, -1.0);

	Camera lCamera(lTopLeftCorner, lImagePlaneWidth, lImagePlaneHeight);

	Vector3d RedColor(1.0, 0.0, 0.0);
	Vector3d BlueColor(0.0, 0.0, 1.0);
	Vector3d GreenColor(0.0, 1.0, 0.0);
	Vector3d YellowColor(1.0, 1.0, 0.4);

	Shapes* lObjlist[3];

	Vector3d lSphere1Center(0.5, 1.0, 0.0);
	//Vector3d lSphere1Center(1.1, 1.1, 0.0);
	Vector3d lSphere2Center(-1.2, 0.9, 0.0);
	Vector3d lSphere3Center(-0.4, 0.0, 0.0);

	ObjType lfullDiffuse(1.0, 0.0, 0.0);
	ObjType lMoreSpecular(0.2, 0.4, 1.0);
	ObjType lNormal;

	lObjlist[0] = new Sphere(lSphere1Center, 0.5, RedColor, lNormal);
	lObjlist[1] = new Sphere(lSphere2Center, 0.3, GreenColor, lfullDiffuse);
	lObjlist[2] = new Sphere(lSphere3Center, 0.7, YellowColor, lMoreSpecular);

	Shapes* lWorld = new ShapesList(lObjlist, 3);

	Vector3d lLightSources[2];
	Vector3d lLightSource1(-1.0, 2.0, -1.0);
	lLightSources[0] = lLightSource1;
	Vector3d lLightSource2(-2.0, -2.0, 1.0);
	lLightSources[1] = lLightSource2;


	Vector3d lEyeSource(-2.0, 0.0, -2.0);

	for (int i = 0; i < liWidth; i++)
	{
		for (int j = 0; j < liHeight; j++)
		{
			//Vector3d lColor(0, 0, 0);

			float lfPixPosX = float(i) / float(liWidth);
			float lfPixPosY = float(j) / float(liHeight);

			Ray lRay = lCamera.getRayP1(lfPixPosX, lfPixPosY, lDirection);

			//for (int k = 0; k < liAvg; k++)
			//{
			//	lColor = lColor + getColorP1(lRay, lWorld);
			//}

			//Vector3d lFColor = lColor / liAvg;

			//Vector3d lFColor = lColor * RedColor;
			Vector3d lFColor = getColorP2(lRay, lWorld, lLightSources, 2, lEyeSource);

			lColorRedValue(i, j) = lFColor[0];
			lColorGreenValue(i, j) = lFColor[1];
			lColorBlueValue(i, j) = lFColor[2];

			// Disable the alpha mask for this pixel
			lAlphaValue(i, j) = 1;
		}

	}


	// Write it in a png image. Note that the alpha channel is reversed to make the white (color = 1) pixels transparent (alhpa = 0)
	write_matrix_to_png(lColorRedValue, lColorGreenValue, lColorBlueValue, lAlphaValue, filename);

}





void part3()
{
	std::cout << "Perspective Projection" << std::endl;

	const std::string filename("part3.png");

	int liWidth = 800;
	int liHeight = 800;
	//int liAvg = 50;

	MatrixXd lColorRedValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorGreenValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorBlueValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lAlphaValue = MatrixXd::Zero(liWidth, liHeight); // Store the alpha mask

	Vector3d lImagePlaneWidth(4.0, 0.0, 0.0);
	Vector3d lImagePlaneHeight(0.0, -4.0, 0.0);
	Vector3d lTopLeftCorner(-2.0, 2.0, 1.0);
	//Vector3d lDirection(0, 0, -1.0);
	Vector3d lOrigin(0.0, 0.4, -1.0);

	Camera lCamera(lTopLeftCorner, lImagePlaneWidth, lImagePlaneHeight);

	Vector3d RedColor(1.0, 0.0, 0.0);
	Vector3d BlueColor(0.0, 0.0, 1.0);
	Vector3d GreenColor(0.0, 1.0, 0.0);
	Vector3d YellowColor(1.0, 1.0, 0.4);

	Shapes* lObjlist[3];

	Vector3d lSphere1Center(0.5, 1.0, 0.0);
	//Vector3d lSphere1Center(1.1, 1.1, 0.0);
	Vector3d lSphere2Center(-1.2, 0.9, 0.0);
	Vector3d lSphere3Center(-0.4, 0.0, 0.0);

	ObjType lfullDiffuse(1.0, 0.0, 0.0);
	ObjType lMoreSpecular(0.2, 0.4, 1.0);
	ObjType lNormal;

	lObjlist[0] = new Sphere(lSphere1Center, 0.5, RedColor, lNormal);
	lObjlist[1] = new Sphere(lSphere2Center, 0.3, GreenColor, lfullDiffuse);
	lObjlist[2] = new Sphere(lSphere3Center, 0.7, YellowColor, lMoreSpecular);

	Shapes* lWorld = new ShapesList(lObjlist, 3);

	Vector3d lLightSources[2];
	Vector3d lLightSource1(-1.0, 2.0, -1.0);
	lLightSources[0] = lLightSource1;
	Vector3d lLightSource2(-2.0, -2.0, 1.0);
	lLightSources[1] = lLightSource2;


	Vector3d lEyeSource(-2.0, 0.0, -2.0);

	for (int i = 0; i < liWidth; i++)
	{
		for (int j = 0; j < liHeight; j++)
		{
			//Vector3d lColor(0, 0, 0);

			float lfPixPosX = float(i) / float(liWidth);
			float lfPixPosY = float(j) / float(liHeight);

			Ray lRay = lCamera.getRayP3(lfPixPosX, lfPixPosY, lOrigin);

			//for (int k = 0; k < liAvg; k++)
			//{
			//	lColor = lColor + getColorP1(lRay, lWorld);
			//}

			//Vector3d lFColor = lColor / liAvg;

			//Vector3d lFColor = lColor * RedColor;
			Vector3d lFColor = getColorP2(lRay, lWorld, lLightSources, 2, lEyeSource);

			lColorRedValue(i, j) = lFColor[0];
			lColorGreenValue(i, j) = lFColor[1];
			lColorBlueValue(i, j) = lFColor[2];

			// Disable the alpha mask for this pixel
			lAlphaValue(i, j) = 1;
		}

	}


	// Write it in a png image. Note that the alpha channel is reversed to make the white (color = 1) pixels transparent (alhpa = 0)
	write_matrix_to_png(lColorRedValue, lColorGreenValue, lColorBlueValue, lAlphaValue, filename);

}

















MatrixXd readMatrixV(const char* pFilename)
{
	int liRows = 0;
	double ldBuff[5];
	int liBuffno = 0;
	int liVSize = 0;
	int liFSize = 0;

	// Read numbers from file into buffer.
	ifstream lInfile;
	lInfile.open(pFilename);

		std::cout  << lInfile.is_open() << std::endl;

	while (!lInfile.eof())
	{
		liBuffno = 0;
		string lLine;
		getline(lInfile, lLine);
		if (liRows == 0)
		{
			liRows++;
			continue;
		}

		if (liRows == 1)
		{
			stringstream stream(lLine);
			while (!stream.eof())
			{
				stream >> ldBuff[liBuffno];
				liBuffno++;
			}
			liVSize = ldBuff[0];
			liFSize = ldBuff[1];
			liRows++;
			//	std::cout << "VSize" << VSize << std::endl;
			//	std::cout << "FSize" << FSize << std::endl;
			continue;
		}
		break;
	}

	lInfile.seekg(0, ios::beg);

	liRows = 0;
	int i = 0;
	MatrixXd result(liVSize, 3);
	while (!lInfile.eof())
	{

		liBuffno = 0;
		string lLine;
		getline(lInfile, lLine);
		//std::cout << lLine << std::endl;

		if (liRows == 0 || liRows == 1 || (liRows > (liVSize + 1)))
		{
			liRows++;
			continue;
		}

		stringstream stream(lLine);
		while (!stream.eof())
		{
			stream >> ldBuff[liBuffno];
			liBuffno++;
		}
		result(i, 0) = ldBuff[0];
		result(i, 1) = ldBuff[1];
		result(i, 2) = ldBuff[2];
		i++;
		liRows++;
	}

	lInfile.close();

	return result;
};


MatrixXi readMatrixF(const char* pFilename)
{
	int liRows = 0;
	double ldBuff[5];
	int liBuffno = 0;
	int liVSize = 0;
	int liFSize = 0;

	// Read numbers from file into buffer.
	ifstream lInfile;
	lInfile.open(pFilename);


	while (!lInfile.eof())
	{
		liBuffno = 0;
		string lLine;
		getline(lInfile, lLine);
		if (liRows == 0)
		{
			liRows++;
			continue;
		}

		if (liRows == 1)
		{
			stringstream stream(lLine);
			while (!stream.eof())
			{
				stream >> ldBuff[liBuffno];
				liBuffno++;
			}
			liVSize = ldBuff[0];
			liFSize = ldBuff[1];
			liRows++;
			std::cout << "VSize" << liVSize << std::endl;
			std::cout << "FSize" << liFSize << std::endl;
			continue;
		}
		break;
	}

	lInfile.seekg(0, ios::beg);

	liRows = 0;
	int i = 0;
	MatrixXi result(liFSize, 4);
	while (!lInfile.eof())
	{

		liBuffno = 0;
		string lLine;
		getline(lInfile, lLine);

		if (liRows < (liVSize + 2) || liRows >(liVSize + 1 + liFSize))
		{
			liRows++;
			continue;
		}

		stringstream stream(lLine);
		while (!stream.eof())
		{
			stream >> ldBuff[liBuffno];
			liBuffno++;
		}
		result(i, 0) = ldBuff[0];
		result(i, 1) = ldBuff[1];
		result(i, 2) = ldBuff[2];
		result(i, 3) = ldBuff[3];
		i++;
		liRows++;
	}

	lInfile.close();

	return result;
};



bool intersectTriangle(Ray pRay, Vector3d pV0, Vector3d pV1, Vector3d pV2, double& pT, double& pU, double& pV)
{
	Vector3d lV1toV0 = pV1 - pV0;
	Vector3d lV2tov0 = pV2 - pV0;

	Vector3d lN = lV1toV0.cross(lV2tov0);
	double lDenominator = lN.dot(lN);

	// ray and plane are parallel
	double lNdotRayDirection = lN.dot(pRay.getDirection());


	//if (fabs(NdotRayDirection) < 0.00000000001) // almost 0 
	//	return false; // they are parallel so they don't intersect ! 

	// compute d parameter using equation 2
	double lD = lN.dot(pV0);

	// compute t (equation 3)
	pT = ( lN.dot(pRay.getOrigin()) + lD) / lNdotRayDirection;

	// check if the triangle is in behind the ray
	if ( pT < 0 ) return false; // the triangle is behind 

	// compute the intersection point using equation 1
	Vector3d lP = pRay.AtPoint(pT);

	// Step 2: inside-outside test
	Vector3d lC; // vector perpendicular to triangle's plane 

	// edge 0
	Vector3d lEdge0 = pV1 - pV0;
	Vector3d lVP0 = lP - pV0;
	lC = lEdge0.cross(lVP0);
	if ( lN.dot(lC) < 0) return false; // P is on the right side 

	// edge 1
	Vector3d lEdge1 = pV2 - pV1;
	Vector3d lVP1 = lP - pV1;
	lC = lEdge1.cross(lVP1);
	if ((pU = lN.dot(lC)) < 0)  return false; // P is on the right side 

	// edge 2
	Vector3d lEdge2 = pV0 - pV2;
	Vector3d lVP2 = lP - pV2;
	lC = lEdge2.cross(lVP2);
	if ((pV = lN.dot(lC)) < 0) return false; // P is on the right side; 

	pU = pU / lDenominator;
	pV = pV / lDenominator;

	return true; // this ray hits the triangle 



}






















void part4()
{
	std::cout << "Ray Tracing Triangle Meshes " << std::endl;

	const std::string filename("part43.png");

	int liWidth = 100;
	int liHeight = 100;
	//int liAvg = 50;

	MatrixXd lColorRedValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorGreenValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lColorBlueValue = MatrixXd::Zero(liWidth, liHeight);
	MatrixXd lAlphaValue = MatrixXd::Zero(liWidth, liHeight); // Store the alpha mask

	Vector3d lImagePlaneWidth(4.0, 0.0, 0.0);
	Vector3d lImagePlaneHeight(0.0, -4.0, 0.0);
	Vector3d lTopLeftCorner(-1.0, 1.0, 1.0);
	//Vector3d lDirection(0, 0, -1.0);
	Vector3d lOrigin(0.0, 0.0, -1.0);

	Camera lCamera(lTopLeftCorner, lImagePlaneWidth, lImagePlaneHeight);

	ifstream fs;

//	fs.open("bumpy_cube.off");
	fs.open("bunny.off");
	if (fs.fail()) {
		std::cout << " bunny.off file not found" << std::endl;
		fs.close();
		return;
	}
	fs.close();


//		MatrixXd lV = readMatrixV("bumpy_cube.off");
//	MatrixXi lF = readMatrixF("bumpy_cube.off");
	MatrixXd lV = readMatrixV("bunny.off");
	MatrixXi lF = readMatrixF("bunny.off");

	
	int liTotalFaces = lF.rows();

	Vector3d lRedColor(1.0, 0.0, 0.0);
	Vector3d lYellowColor(1.0, 1.0, 0.0);
	Vector3d lGreenColor(0.0, 1.0, 0.0);
	Vector3d lBlueColor(0.0, 0.0, 1.0);




	for (int i = 0; i < liWidth; i++)
	{
		for (int j = 0; j < liHeight; j++)
		{


			std::cout << " i " << i << " j " << j << std::endl;

			float lfPixPosX = float( (i - (liWidth/2.0) )/ float(liWidth) );
			float lfPixPosY = float( (j - (liHeight/2.0) ) / float(liHeight) );

			Ray lRay = lCamera.getRayP3(lfPixPosX, lfPixPosY, lOrigin);


			double lT,lPU,lPV,lLU = 1.0, lLV = 1.0 ,lLT=getMaxFloatNumber();
			int liIntersect = 0;

			for (int f = 0; f < liTotalFaces; f++)
			{
				Vector3d v0(lV(lF(f, 1), 0), lV(lF(f, 1), 1), lV(lF(f, 1), 2));
				Vector3d v1(lV(lF(f, 2), 0), lV(lF(f, 2), 1), lV(lF(f, 2), 2));
				Vector3d v2(lV(lF(f, 3), 0), lV(lF(f, 3), 1), lV(lF(f, 3), 2));

				if ( intersectTriangle(lRay, v0, v1, v2, lT, lPU, lPV) )
				{
					if ( lT < lLT && lT > 0)
					{
						lLT = lT;
						lLU = lPU;
						lLV = lPV;
						liIntersect = 1;
					}
				}

			}

			Vector3d lFColor;

			if (liIntersect == 1)
			{
				lLU = clampNumberBetween0and1(lLU);
				lLV = clampNumberBetween0and1(lLV);
				lFColor = (lLU * lRedColor) + (lLV * lYellowColor) + ((1 - lLU - lLV) * lGreenColor);
			}
			else
			{
				lFColor = lBlueColor;
			}


			lColorRedValue(i, j) = lFColor[0];
			lColorGreenValue(i, j) = lFColor[1];
			lColorBlueValue(i, j) = lFColor[2];

			// Disable the alpha mask for this pixel
			lAlphaValue(i, j) = 1;
		}

	}

	

	// Write it in a png image. Note that the alpha channel is reversed to make the white (color = 1) pixels transparent (alhpa = 0)
	write_matrix_to_png(lColorRedValue, lColorGreenValue, lColorBlueValue, lAlphaValue, filename);


}


int main()
{
   // part1();
   // part2();
//	part3();
	part4();

    return 0;
}
