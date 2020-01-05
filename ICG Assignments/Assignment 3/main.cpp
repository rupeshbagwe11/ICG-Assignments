// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif



#include <iostream>
#include <string>
#include <cstdlib>
#include<cmath>
#include <limits>
#include <fstream>
#include <list> 
#include <vector> 
#include <iterator> 
#include <algorithm>


#include <Eigen/Core>
#include <Eigen/Geometry>


using namespace std;
using namespace Eigen;

// Timer
#include <chrono>

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_C;
VertexBufferObject VBO_NF;
VertexBufferObject VBO_NP;

// Contains the vertex positions
Eigen::MatrixXf V(3, 1);
Eigen::MatrixXf C(3, 1);
Eigen::MatrixXf NF(3, 1);
Eigen::MatrixXf NP(3, 1);


//sadsads

int liBuTotalVertex = 0;
int liBcTotalVertex = 0;
int liCubeTotalVertex = 36;

int miCurrentMode = 1; //i = 1, o = 2, c = 3, b = 4
int miProjectionMode = 1; //z = 1, x = 2
int miCameraMode = 1; //normal = 1, circular = 2
int miNoOfObjects = 0;

int miSelectedObject = -1;
int mbIsPressed = 0;

int miWidth, miHeight;
double tempno = 0;

double mdSelectedObjectX = 0.0;
double mdSelectedObjectY = 0.0;

vector<int> mlistObjType, mlistObjDispType, mlistObjR, mlistObjS, mlistObjDeleted;
vector<float>  mlistObjTX, mlistObjTY;

double mdCameraX = 0, mdCameraY = 0, mdCameraZ = 1;
int miXPress = 0, miYPress = 0, miZPress = 0;
double PI = 3.14159;

Vector3d mLPinkColor(1.0f, 0.2f, 0.6f);
Vector3d mLBlueColor(0.2f, 0.6f, 1.0f);
Vector3d mWhiteColor(1.0f, 1.0f, 1.0f);
Vector3d mDGreenColor(0.0f, 0.8f, 0.8f);
Vector3d mOrangeColor(1.0f, 0.5f, 0.0f);
Vector3d mBlackColor(0.0f, 0.0f, 0.0f);
Vector3d mYellowColor(1.0f, 1.0f, 0.2f);
Vector3d mRedColor(1.0f, 0.0f, 0.0f);
Vector3d mBlueColor(0.0f, 0.0f, 1.0f);
Vector3d mGreenColor(0.0f, 1.0f, 0.0f);



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
			//		std::cout << "VSize" << liVSize << std::endl;
			//		std::cout << "FSize" << liFSize << std::endl;
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



MatrixXd readMatrixV(const char* pFilename, double* lVxmin, double* lVxmax, double* lVymin, double* lVymax, double* lVzmin, double* lVzmax)
{
	int liRows = 0;
	double ldBuff[5];
	int liBuffno = 0;
	int liVSize = 0;
	int liFSize = 0;
	double lxmin = 99999999999999;
	double lxmax = -99999999999999;

	double lymin = 99999999999999;
	double lymax = -99999999999999;

	double lzmin = 99999999999999;
	double lzmax = -99999999999999;

	// Read numbers from file into buffer.
	ifstream lInfile;
	lInfile.open(pFilename);

	//std::cout << lInfile.is_open() << std::endl;

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

		if (ldBuff[0] < lxmin)
		{
			lxmin = ldBuff[0];
		}
		if (ldBuff[1] < lymin)
		{
			lymin = ldBuff[1];
		}
		if (ldBuff[2] < lzmin)
		{
			lzmin = ldBuff[2];
		}

		if (ldBuff[0] > lxmax)
		{
			lxmax = ldBuff[0];
		}
		if (ldBuff[1] > lymax)
		{
			lymax = ldBuff[1];
		}
		if (ldBuff[2] > lzmax)
		{
			lzmax = ldBuff[2];
		}

		result(i, 0) = ldBuff[0];
		result(i, 1) = ldBuff[1];
		result(i, 2) = ldBuff[2];
		i++;
		liRows++;
	}

	lInfile.close();
	*lVxmin = lxmin;
	*lVxmax = lxmax;
	*lVymin = lymin;
	*lVymax = lymax;
	*lVzmin = lzmin;
	*lVzmax = lzmax;

	return result;
};


double convertToRange(double omin, double omax, double nmin, double nmax, double no)
{
	double or = (omax - omin);
	double nr = (nmax - nmin);
	double nv = (((no - omin) * nr) / or) + nmin;
	return nv;
}

//rteun the start index for the next object if objectno = -1 else the given objectno
int getStartIndex(int pObjectNo)
{
	int liSIndex = 0;
	int lisize = 0;
	if (pObjectNo == -1)
	{
		lisize = mlistObjType.size();
	}
	else
	{
		lisize = pObjectNo - 1;
	}
	for (int i = 0; i < lisize; i++)
	{
		if (mlistObjType[i] == 1)
		{
			liSIndex = liSIndex + liCubeTotalVertex;
		}
		else if (mlistObjType[i] == 2)
		{
			liSIndex = liSIndex + liBuTotalVertex;
		}
		else if (mlistObjType[i] == 3)
		{
			liSIndex = liSIndex + liBcTotalVertex;
		}
	}
	return liSIndex;
}



void addBunny()
{
	ifstream fs;

	fs.open("bunny.off");
	if (fs.fail()) {
		std::cout << " bunny.off file not found" << std::endl;
		fs.close();
		return;
	}
	fs.close();

	double lVxmin = 0, lVxmax = 0, lVymin = 0, lVymax = 0, lVzmin = 0, lVzmax = 0;
	//std::cout << lVmin << " " << lVmax << "\n";
	MatrixXd lBuV = readMatrixV("bunny.off", &lVxmin, &lVxmax, &lVymin, &lVymax, &lVzmin, &lVzmax);
	MatrixXi lBuF = readMatrixF("bunny.off");

	int liBuTotalFaces = lBuF.rows();
	liBuTotalVertex = liBuTotalFaces * 3;
	//std::cout << lVmin << " " << lVmax << "\n";

	V.conservativeResize(V.rows(), V.cols() + liBuTotalVertex);
	C.conservativeResize(C.rows(), C.cols() + liBuTotalVertex);
	NF.conservativeResize(NF.rows(), NF.cols() + liBuTotalVertex);
	NP.conservativeResize(NP.rows(), NP.cols() + liBuTotalVertex);
	//std::cout << lBuV << "\n";
	//std::cout << lBuF << "\n";

	int sIndex = getStartIndex(-1);
	for (int f = 0; f < liBuTotalFaces; f++)
	{
		V.col(sIndex + ((f * 3) + 0)) << convertToRange(lVxmin, lVxmax, -0.5, 0.5, lBuV(lBuF(f, 1), 0)), convertToRange(lVymin, lVymax, -0.5, 0.5, lBuV(lBuF(f, 1), 1)), convertToRange(lVzmin, lVzmax, -0.5, 0.5, lBuV(lBuF(f, 1), 2));
		V.col(sIndex + ((f * 3) + 1)) << convertToRange(lVxmin, lVxmax, -0.5, 0.5, lBuV(lBuF(f, 2), 0)), convertToRange(lVymin, lVymax, -0.5, 0.5, lBuV(lBuF(f, 2), 1)), convertToRange(lVzmin, lVzmax, -0.5, 0.5, lBuV(lBuF(f, 2), 2));
		V.col(sIndex + ((f * 3) + 2)) << convertToRange(lVxmin, lVxmax, -0.5, 0.5, lBuV(lBuF(f, 3), 0)), convertToRange(lVymin, lVymax, -0.5, 0.5, lBuV(lBuF(f, 3), 1)), convertToRange(lVzmin, lVzmax, -0.5, 0.5, lBuV(lBuF(f, 3), 2));

		C.col(sIndex + ((f * 3) + 0)) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
		C.col(sIndex + ((f * 3) + 1)) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
		C.col(sIndex + ((f * 3) + 2)) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);

	}

	for (int j = 0; j < 3000; j = j + 3)
	{
		Vector3f lV1 = V.col(sIndex + j + 0);
		Vector3f lV2 = V.col(sIndex + j + 1);
		Vector3f lV3 = V.col(sIndex + j + 2);
		Vector3f lU = lV2 - lV1;
		Vector3f lV = lV3 - lV1;
		Vector3f lNormal;
		lNormal << (lU(1) * lV(2)) - (lU(2) * lV(1)), ((lU(2) * lV(0)) - (lU(0) * lV(2))), ((lU(0) * lV(1)) - (lU(1) * lV(0)));
		float ldist = sqrt((lNormal(0) * lNormal(0)) + (lNormal(1) * lNormal(1)) + (lNormal(2) * lNormal(2)));

		NF.col(sIndex + j + 0) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
		NF.col(sIndex + j + 1) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
		NF.col(sIndex + j + 2) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
	}

	Eigen::MatrixXf NV(4, 502);
	for (int f = 0; f < liBuTotalFaces; f++)
	{
		Vector4f lTempV;
		Vector4f lfV;

		float lVertexIndex = lBuF(f, 1);
		Vector3f lTempN = NF.col(sIndex + (f * 3) + 0);
		std::cout << std::fixed;
		//	std::cout << NV.col(0)(2) << "\n";
		if (NV.col(lVertexIndex)(3) < -100)
		{
			//		std::cout << NV.col(0)(3) << "less \n";
			lTempV << NV.col(lVertexIndex)(0), NV.col(lVertexIndex)(1), NV.col(lVertexIndex)(2), NV.col(lVertexIndex)(3);
			lfV << lTempN(0), lTempN(1), lTempN(2), 1;
		}
		else
		{
			//		std::cout << NV.col(0)(3) << "more \n";
			lTempV << NV.col(lVertexIndex)(0), NV.col(lVertexIndex)(1), NV.col(lVertexIndex)(2), NV.col(lVertexIndex)(3);
			lfV << lTempV(0) + lTempN(0), lTempV(1) + lTempN(1), lTempV(2) + lTempN(2), lTempV(3) + 1;
		}
		//std::cout << lTempV << "\n";
		NV.col(lVertexIndex) = lfV;

		float lVertexIndex1 = lBuF(f, 2);
		Vector3f lTempN1 = NF.col(sIndex + (f * 3) + 1);
		Vector4f lTempV1;
		Vector4f lfV1;

		if (NV.col(lVertexIndex1)(3) < -100)
		{
			lTempV1 << NV.col(lVertexIndex1)(0), NV.col(lVertexIndex1)(1), NV.col(lVertexIndex1)(2), NV.col(lVertexIndex1)(3);
			lfV1 << lTempN1(0), lTempN1(1), lTempN1(2), 1;
		}
		else
		{
			lTempV1 << NV.col(lVertexIndex1)(0), NV.col(lVertexIndex1)(1), NV.col(lVertexIndex1)(2), NV.col(lVertexIndex1)(3);
			lfV1 << lTempV1(0) + lTempN1(0), lTempV1(1) + lTempN1(1), lTempV1(2) + lTempN1(2), lTempV1(3) + 1;
		}
		NV.col(lVertexIndex1) = lfV;


		float lVertexIndex2 = lBuF(f, 3);
		Vector3f lTempN2 = NF.col(sIndex + (f * 3) + 2);
		Vector4f lTempV2;
		Vector4f lfV2;

		if (NV.col(lVertexIndex2)(3) < -100)
		{
			lTempV2 << NV.col(lVertexIndex2)(0), NV.col(lVertexIndex2)(1), NV.col(lVertexIndex2)(2), NV.col(lVertexIndex2)(3);
			lfV2 << lTempN2(0), lTempN2(1), lTempN2(2), 1;
		}
		else
		{
			lTempV2 << NV.col(lVertexIndex2)(0), NV.col(lVertexIndex2)(1), NV.col(lVertexIndex2)(2), NV.col(lVertexIndex2)(3);
			lfV2 << lTempV2(0) + lTempN2(0), lTempV2(1) + lTempN2(1), lTempV2(2) + lTempN2(2), lTempV2(3) + 1;
		}
		NV.col(lVertexIndex2) = lfV;
	}

	for (int f = 0; f < liBuTotalFaces; f++)
	{
		float lVertexIndex = lBuF(f, 1);
		Vector3f lNormal;
		lNormal << NV.col(lVertexIndex)(0) / NV.col(lVertexIndex)(3), NV.col(lVertexIndex)(1) / NV.col(lVertexIndex)(3), NV.col(lVertexIndex)(2) / NV.col(lVertexIndex)(3);
		float ldist = sqrt((lNormal(0) * lNormal(0)) + (lNormal(1) * lNormal(1)) + (lNormal(2) * lNormal(2)));
		NP.col(sIndex + (f * 3 + 0)) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;

		float lVertexIndex1 = lBuF(f, 2);
		Vector3f lNormal1;
		lNormal1 << NV.col(lVertexIndex1)(0) / NV.col(lVertexIndex1)(3), NV.col(lVertexIndex1)(1) / NV.col(lVertexIndex1)(3), NV.col(lVertexIndex1)(2) / NV.col(lVertexIndex1)(3);
		float ldist1 = sqrt((lNormal1(0) * lNormal1(0)) + (lNormal1(1) * lNormal1(1)) + (lNormal1(2) * lNormal1(2)));
		NP.col(sIndex + (f * 3 + 1)) << lNormal1(0) / ldist1, lNormal1(1) / ldist1, lNormal1(2) / ldist1;

		float lVertexIndex2 = lBuF(f, 3);
		Vector3f lNormal2;
		lNormal2 << NV.col(lVertexIndex2)(0) / NV.col(lVertexIndex2)(3), NV.col(lVertexIndex2)(1) / NV.col(lVertexIndex2)(3), NV.col(lVertexIndex2)(2) / NV.col(lVertexIndex2)(3);
		float ldist2 = sqrt((lNormal2(0) * lNormal2(0)) + (lNormal2(1) * lNormal2(1)) + (lNormal2(2) * lNormal2(2)));
		NP.col(sIndex + (f * 3 + 2)) << lNormal2(0) / ldist2, lNormal2(1) / ldist2, lNormal2(2) / ldist2;

	}


	miNoOfObjects++;
	mlistObjType.push_back(2);
	mlistObjDispType.push_back(3);

	mlistObjTX.push_back(0.0f);
	mlistObjTY.push_back(0.0f);
	mlistObjR.push_back(0);
	mlistObjS.push_back(0);
	mlistObjDeleted.push_back(0);

	//std::cout << V << "\n";
	VBO.update(V);
	VBO_C.update(C);
	VBO_NF.update(NF);
	VBO_NP.update(NP);
}

void addBumpyCube()
{
	ifstream fs;

	fs.open("bumpy_cube.off");
	if (fs.fail()) {
		std::cout << " bumpy_cube.off file not found" << std::endl;
		fs.close();
		return;
	}
	fs.close();

	double lVxmin = 0, lVxmax = 0, lVymin = 0, lVymax = 0, lVzmin = 0, lVzmax = 0;
	//std::cout << lVmin << " " << lVmax << "\n";
	MatrixXd lBcV = readMatrixV("bumpy_cube.off", &lVxmin, &lVxmax, &lVymin, &lVymax, &lVzmin, &lVzmax);
	MatrixXi lBcF = readMatrixF("bumpy_cube.off");

	int liBcTotalFaces = lBcF.rows();
	liBcTotalVertex = liBcTotalFaces * 3;
	//std::cout << lVmin << " " << lVmax << "\n";

	V.conservativeResize(V.rows(), V.cols() + liBcTotalVertex);
	C.conservativeResize(C.rows(), C.cols() + liBcTotalVertex);
	NF.conservativeResize(NF.rows(), NF.cols() + liBcTotalVertex);
	NP.conservativeResize(NP.rows(), NP.cols() + liBcTotalVertex);
	//std::cout << lBuV << "\n";
	//std::cout << lBuF << "\n";

	int sIndex = getStartIndex(-1);
	for (int f = 0; f < liBcTotalFaces; f++)
	{
		V.col(sIndex + ((f * 3) + 0)) << convertToRange(lVxmin, lVxmax, -0.5, 0.5, lBcV(lBcF(f, 1), 0)), convertToRange(lVymin, lVymax, -0.5, 0.5, lBcV(lBcF(f, 1), 1)), convertToRange(lVzmin, lVzmax, -0.5, 0.5, lBcV(lBcF(f, 1), 2));
		V.col(sIndex + ((f * 3) + 1)) << convertToRange(lVxmin, lVxmax, -0.5, 0.5, lBcV(lBcF(f, 2), 0)), convertToRange(lVymin, lVymax, -0.5, 0.5, lBcV(lBcF(f, 2), 1)), convertToRange(lVzmin, lVzmax, -0.5, 0.5, lBcV(lBcF(f, 2), 2));
		V.col(sIndex + ((f * 3) + 2)) << convertToRange(lVxmin, lVxmax, -0.5, 0.5, lBcV(lBcF(f, 3), 0)), convertToRange(lVymin, lVymax, -0.5, 0.5, lBcV(lBcF(f, 3), 1)), convertToRange(lVzmin, lVzmax, -0.5, 0.5, lBcV(lBcF(f, 3), 2));

		C.col(sIndex + ((f * 3) + 0)) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
		C.col(sIndex + ((f * 3) + 1)) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
		C.col(sIndex + ((f * 3) + 2)) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
	}

	for (int j = 0; j < 3000; j = j + 3)
	{
		Vector3f lV1 = V.col(sIndex + j + 0);
		Vector3f lV2 = V.col(sIndex + j + 1);
		Vector3f lV3 = V.col(sIndex + j + 2);
		Vector3f lU = lV2 - lV1;
		Vector3f lV = lV3 - lV1;
		Vector3f lNormal;
		lNormal << (lU(1) * lV(2)) - (lU(2) * lV(1)), ((lU(2) * lV(0)) - (lU(0) * lV(2))), ((lU(0) * lV(1)) - (lU(1) * lV(0)));
		float ldist = sqrt((lNormal(0) * lNormal(0)) + (lNormal(1) * lNormal(1)) + (lNormal(2) * lNormal(2)));

		NF.col(sIndex + j + 0) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
		NF.col(sIndex + j + 1) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
		NF.col(sIndex + j + 2) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
	}

	Eigen::MatrixXf NV(4, 502);
	for (int f = 0; f < liBcTotalFaces; f++)
	{
		Vector4f lTempV;
		Vector4f lfV;

		float lVertexIndex = lBcF(f, 1);
		Vector3f lTempN = NF.col(sIndex + (f * 3) + 0);
		//std::cout << std::fixed;
		//	std::cout << NV.col(0)(2) << "\n";
		if (NV.col(lVertexIndex)(3) < -100)
		{
			//		std::cout << NV.col(0)(3) << "less \n";
			lTempV << NV.col(lVertexIndex)(0), NV.col(lVertexIndex)(1), NV.col(lVertexIndex)(2), NV.col(lVertexIndex)(3);
			lfV << lTempN(0), lTempN(1), lTempN(2), 1;
		}
		else
		{
			//		std::cout << NV.col(0)(3) << "more \n";
			lTempV << NV.col(lVertexIndex)(0), NV.col(lVertexIndex)(1), NV.col(lVertexIndex)(2), NV.col(lVertexIndex)(3);
			lfV << lTempV(0) + lTempN(0), lTempV(1) + lTempN(1), lTempV(2) + lTempN(2), lTempV(3) + 1;
		}
		//std::cout << lTempV << "\n";
		NV.col(lVertexIndex) = lfV;

		float lVertexIndex1 = lBcF(f, 2);
		Vector3f lTempN1 = NF.col(sIndex + (f * 3) + 1);
		Vector4f lTempV1;
		Vector4f lfV1;

		if (NV.col(lVertexIndex1)(3) < -100)
		{
			lTempV1 << NV.col(lVertexIndex1)(0), NV.col(lVertexIndex1)(1), NV.col(lVertexIndex1)(2), NV.col(lVertexIndex1)(3);
			lfV1 << lTempN1(0), lTempN1(1), lTempN1(2), 1;
		}
		else
		{
			lTempV1 << NV.col(lVertexIndex1)(0), NV.col(lVertexIndex1)(1), NV.col(lVertexIndex1)(2), NV.col(lVertexIndex1)(3);
			lfV1 << lTempV1(0) + lTempN1(0), lTempV1(1) + lTempN1(1), lTempV1(2) + lTempN1(2), lTempV1(3) + 1;
		}
		NV.col(lVertexIndex1) = lfV;


		float lVertexIndex2 = lBcF(f, 3);
		Vector3f lTempN2 = NF.col(sIndex + (f * 3) + 2);
		Vector4f lTempV2;
		Vector4f lfV2;

		if (NV.col(lVertexIndex2)(3) < -100)
		{
			lTempV2 << NV.col(lVertexIndex2)(0), NV.col(lVertexIndex2)(1), NV.col(lVertexIndex2)(2), NV.col(lVertexIndex2)(3);
			lfV2 << lTempN2(0), lTempN2(1), lTempN2(2), 1;
		}
		else
		{
			lTempV2 << NV.col(lVertexIndex2)(0), NV.col(lVertexIndex2)(1), NV.col(lVertexIndex2)(2), NV.col(lVertexIndex2)(3);
			lfV2 << lTempV2(0) + lTempN2(0), lTempV2(1) + lTempN2(1), lTempV2(2) + lTempN2(2), lTempV2(3) + 1;
		}
		NV.col(lVertexIndex2) = lfV;
	}

	for (int f = 0; f < liBcTotalFaces; f++)
	{
		float lVertexIndex = lBcF(f, 1);
		Vector3f lNormal;
		lNormal << NV.col(lVertexIndex)(0) / NV.col(lVertexIndex)(3), NV.col(lVertexIndex)(1) / NV.col(lVertexIndex)(3), NV.col(lVertexIndex)(2) / NV.col(lVertexIndex)(3);
		float ldist = sqrt((lNormal(0) * lNormal(0)) + (lNormal(1) * lNormal(1)) + (lNormal(2) * lNormal(2)));
		NP.col(sIndex + (f * 3 + 0)) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;

		float lVertexIndex1 = lBcF(f, 2);
		Vector3f lNormal1;
		lNormal1 << NV.col(lVertexIndex1)(0) / NV.col(lVertexIndex1)(3), NV.col(lVertexIndex1)(1) / NV.col(lVertexIndex1)(3), NV.col(lVertexIndex1)(2) / NV.col(lVertexIndex1)(3);
		float ldist1 = sqrt((lNormal1(0) * lNormal1(0)) + (lNormal1(1) * lNormal1(1)) + (lNormal1(2) * lNormal1(2)));
		NP.col(sIndex + (f * 3 + 1)) << lNormal1(0) / ldist1, lNormal1(1) / ldist1, lNormal1(2) / ldist1;

		float lVertexIndex2 = lBcF(f, 3);
		Vector3f lNormal2;
		lNormal2 << NV.col(lVertexIndex2)(0) / NV.col(lVertexIndex2)(3), NV.col(lVertexIndex2)(1) / NV.col(lVertexIndex2)(3), NV.col(lVertexIndex2)(2) / NV.col(lVertexIndex2)(3);
		float ldist2 = sqrt((lNormal2(0) * lNormal2(0)) + (lNormal2(1) * lNormal2(1)) + (lNormal2(2) * lNormal2(2)));
		NP.col(sIndex + (f * 3 + 2)) << lNormal2(0) / ldist2, lNormal2(1) / ldist2, lNormal2(2) / ldist2;

	}


	miNoOfObjects++;
	mlistObjType.push_back(3);
	mlistObjDispType.push_back(3);

	mlistObjTX.push_back(0.0f);
	mlistObjTY.push_back(0.0f);
	mlistObjR.push_back(0);
	mlistObjS.push_back(0);
	mlistObjDeleted.push_back(0);

	//std::cout << V << "\n";
	VBO.update(V);
	VBO_C.update(C);
	VBO_NF.update(NF);
	VBO_NP.update(NP);
}



void addCube()
{
	V.conservativeResize(V.rows(), V.cols() + liCubeTotalVertex);
	C.conservativeResize(C.rows(), C.cols() + liCubeTotalVertex);
	NF.conservativeResize(NF.rows(), NF.cols() + liCubeTotalVertex);
	NP.conservativeResize(NP.rows(), NP.cols() + liCubeTotalVertex);

	int sIndex = getStartIndex(-1);

	V.col(sIndex + 0) << -0.5, 0.5, 0.5;//1 -++
	V.col(sIndex + 1) << -0.5, -0.5, 0.5;//2 --+
	V.col(sIndex + 2) << 0.5, -0.5, 0.5;//3 +-+

	V.col(sIndex + 3) << -0.5, 0.5, 0.5;//1 -++
	V.col(sIndex + 4) << 0.5, 0.5, 0.5;//4 +++
	V.col(sIndex + 5) << 0.5, -0.5, 0.5;//3 +-+

	V.col(sIndex + 6) << 0.5, 0.5, 0.5;//4 +++
	V.col(sIndex + 7) << 0.5, -0.5, 0.5;//3 +-+
	V.col(sIndex + 8) << 0.5, -0.5, -0.5;//5 +--

	V.col(sIndex + 9) << 0.5, 0.5, 0.5;//4 +++
	V.col(sIndex + 10) << 0.5, 0.5, -0.5;//6 ++-
	V.col(sIndex + 11) << 0.5, -0.5, -0.5;//5 +--

	V.col(sIndex + 12) << -0.5, 0.5, -0.5;//7 -+-
	V.col(sIndex + 13) << -0.5, -0.5, -0.5;//8 ---
	V.col(sIndex + 14) << 0.5, -0.5, -0.5;//5 +--

	V.col(sIndex + 15) << -0.5, 0.5, -0.5;//7 -+-
	V.col(sIndex + 16) << 0.5, 0.5, -0.5;//6 ++-
	V.col(sIndex + 17) << 0.5, -0.5, -0.5;//5 +--

	V.col(sIndex + 18) << -0.5, 0.5, 0.5;//1 -++
	V.col(sIndex + 19) << -0.5, -0.5, 0.5;//2 --+
	V.col(sIndex + 20) << -0.5, -0.5, -0.5;//8 ---

	V.col(sIndex + 21) << -0.5, 0.5, 0.5;//1 -++
	V.col(sIndex + 22) << -0.5, 0.5, -0.5;//7 -+-
	V.col(sIndex + 23) << -0.5, -0.5, -0.5;//8 ---

	V.col(sIndex + 24) << -0.5, 0.5, -0.5;//7 -+-
	V.col(sIndex + 25) << -0.5, 0.5, 0.5;//1 -++
	V.col(sIndex + 26) << 0.5, 0.5, 0.5;//4 +++

	V.col(sIndex + 27) << -0.5, 0.5, -0.5;//7 -+-
	V.col(sIndex + 28) << 0.5, 0.5, -0.5;//6 ++-
	V.col(sIndex + 29) << 0.5, 0.5, 0.5;//4 +++

	V.col(sIndex + 30) << -0.5, -0.5, -0.5;//8 ---
	V.col(sIndex + 31) << -0.5, -0.5, 0.5;//2 --+
	V.col(sIndex + 32) << 0.5, -0.5, 0.5;//3 +-+

	V.col(sIndex + 33) << -0.5, -0.5, -0.5;//8 ---
	V.col(sIndex + 34) << 0.5, -0.5, -0.5;//5 +--
	V.col(sIndex + 35) << 0.5, -0.5, 0.5;//3 +-+

	for (int j = 0; j < 36; j++)
	{
		C.col(sIndex + j) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
	}

	for (int j = 0; j < 36; j = j + 3)
	{
		Vector3f lV1 = V.col(sIndex + j + 0);
		Vector3f lV2 = V.col(sIndex + j + 1);
		Vector3f lV3 = V.col(sIndex + j + 2);
		Vector3f lU = lV2 - lV1;
		Vector3f lV = lV3 - lV1;
		Vector3f lNormal;
		lNormal << (lU(1) * lV(2)) - (lU(2) * lV(1)), ((lU(2) * lV(0)) - (lU(0) * lV(2))), ((lU(0) * lV(1)) - (lU(1) * lV(0)));
		float ldist = sqrt((lNormal(0) * lNormal(0)) + (lNormal(1) * lNormal(1)) + (lNormal(2) * lNormal(2)));

		NF.col(sIndex + j + 0) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
		NF.col(sIndex + j + 1) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
		NF.col(sIndex + j + 2) << lNormal(0) / ldist, lNormal(1) / ldist, lNormal(2) / ldist;
	}

	NP.col(sIndex + 0) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 25)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 25)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 25)(2) / 5;
	NP.col(sIndex + 1) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 30)(0) / 3, NF.col(sIndex + 1)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 30)(1) / 3, NF.col(sIndex + 1)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 30)(2) / 3;
	NP.col(sIndex + 2) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 3) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 25)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 25)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 25)(2) / 5;
	NP.col(sIndex + 4) << NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 5) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 6) << NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 7) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 8) << NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 9) << NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 10) << NF.col(sIndex + 9)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 27)(0) / 3, NF.col(sIndex + 9)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 27)(1) / 3, NF.col(sIndex + 9)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 27)(2) / 3;
	NP.col(sIndex + 11) << NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 12) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 13) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 14) << NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 15) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 16) << NF.col(sIndex + 9)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 27)(0) / 3, NF.col(sIndex + 9)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 27)(1) / 3, NF.col(sIndex + 9)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 27)(2) / 3;
	NP.col(sIndex + 17) << NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 18) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 25)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 25)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 25)(2) / 5;
	NP.col(sIndex + 19) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 20) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 21) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 25)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 25)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 25)(2) / 5;
	NP.col(sIndex + 22) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 23) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 24) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 25) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 25)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 25)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 25)(2) / 5;
	NP.col(sIndex + 26) << NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;

	NP.col(sIndex + 27) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;
	NP.col(sIndex + 28) << NF.col(sIndex + 9)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 27)(0) / 3, NF.col(sIndex + 9)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 27)(1) / 3, NF.col(sIndex + 9)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 27)(2) / 3;
	NP.col(sIndex + 29) << NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 24)(0) + NF.col(sIndex + 27)(0) / 5, NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 24)(1) + NF.col(sIndex + 27)(1) / 5, NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 24)(2) + NF.col(sIndex + 27)(2) / 5;

	NP.col(sIndex + 30) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 31) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 32) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;

	NP.col(sIndex + 33) << NF.col(sIndex + 12)(0) + NF.col(sIndex + 18)(0) + NF.col(sIndex + 21)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 12)(1) + NF.col(sIndex + 18)(1) + NF.col(sIndex + 21)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 12)(2) + NF.col(sIndex + 18)(2) + NF.col(sIndex + 21)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 34) << NF.col(sIndex + 6)(0) + NF.col(sIndex + 9)(0) + NF.col(sIndex + 12)(0) + NF.col(sIndex + 15)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 6)(1) + NF.col(sIndex + 9)(1) + NF.col(sIndex + 12)(1) + NF.col(sIndex + 15)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 6)(2) + NF.col(sIndex + 9)(2) + NF.col(sIndex + 12)(2) + NF.col(sIndex + 15)(2) + NF.col(sIndex + 33)(2) / 5;
	NP.col(sIndex + 35) << NF.col(sIndex + 1)(0) + NF.col(sIndex + 3)(0) + NF.col(sIndex + 6)(0) + NF.col(sIndex + 30)(0) + NF.col(sIndex + 33)(0) / 5, NF.col(sIndex + 1)(1) + NF.col(sIndex + 3)(1) + NF.col(sIndex + 6)(1) + NF.col(sIndex + 30)(1) + NF.col(sIndex + 33)(1) / 5, NF.col(sIndex + 1)(2) + NF.col(sIndex + 3)(2) + NF.col(sIndex + 6)(2) + NF.col(sIndex + 30)(2) + NF.col(sIndex + 33)(2) / 5;



	miNoOfObjects++;
	mlistObjType.push_back(1);
	mlistObjDispType.push_back(3);

	mlistObjTX.push_back(0.0f);
	mlistObjTY.push_back(0.0f);
	mlistObjR.push_back(0);
	mlistObjS.push_back(0);
	mlistObjDeleted.push_back(0);

	VBO.update(V);
	VBO_C.update(C);
	VBO_NF.update(NF);
	VBO_NP.update(NP);
}

static double area(double pV1x, double pV1y, double pV2x, double pV2y, double pV3x, double pV3y)
{
	return abs((pV1x * (pV2y - pV3y) + pV2x * (pV3y - pV1y) + pV3x * (pV1y - pV2y)) / 2.0);
}

static boolean isInsideTriangle(double pTV1x, double pTV1y, double pTV2x, double pTV2y, double pTV3x, double pTV3y, double pPx, double pPy)
{
	/* Area of ABC */
	double ldAT = area(pTV1x, pTV1y, pTV2x, pTV2y, pTV3x, pTV3y);

	/* Area of PBC */
	double ldAST1 = area(pPx, pPy, pTV2x, pTV2y, pTV3x, pTV3y);

	/* Area of PAC */
	double ldAST2 = area(pTV1x, pTV1y, pPx, pPy, pTV3x, pTV3y);

	/* Area of PAB */
	double ldAST3 = area(pTV1x, pTV1y, pTV2x, pTV2y, pPx, pPy);

	//	std::cout << std::fixed;
	//    std:cout << " ld " << ldAT << " sum" << (ldAST1 + ldAST2 + ldAST3) << " " << (ldAT <= ldAST1 + ldAST2 + ldAST3 + 0.00002) << " " << (ldAT >= ldAST1 + ldAST2 + ldAST3 - 0.00002) <<" \n";
		/* A1+A2+A3 = A with 0.002 calculateError */
	return (ldAT <= ldAST1 + ldAST2 + ldAST3 + 0.00002) && (ldAT >= ldAST1 + ldAST2 + ldAST3 - 0.00002);
}


template<class T>
Eigen::Matrix<T, 4, 4> getPerspectiveMatrix(double pFOV, double pAspect, double pNear, double pFar)
{
	typedef Eigen::Matrix<T, 4, 4> Matrix4;

	double ldRadf = (pFOV * PI) / 180;

	double lTanHalfFOV = tan(ldRadf / 2.0);
	Matrix4 lRes = Matrix4::Zero();
	lRes(0, 0) = 1.0 / (pAspect * lTanHalfFOV);
	lRes(1, 1) = 1.0 / (lTanHalfFOV);
	lRes(2, 2) = -(pFar + pNear) / (pFar - pNear);
	lRes(3, 2) = -1.0;
	lRes(2, 3) = -(2.0 * pFar * pNear) / (pFar - pNear);
	return lRes;
}

template<class T>
Eigen::Matrix<T, 4, 4> getLookAtMatrix
(
	Eigen::Matrix<T, 3, 1> const& pEye,
	Eigen::Matrix<T, 3, 1> const& pCenter,
	Eigen::Matrix<T, 3, 1> const& pUp
)
{
	typedef Eigen::Matrix<T, 4, 4> Matrix4;
	typedef Eigen::Matrix<T, 3, 1> Vector3;

	Vector3 lF = (pCenter - pEye).normalized();
	Vector3 lU = pUp.normalized();
	Vector3 lS = lF.cross(lU).normalized();
	lU = lS.cross(lF);

	Matrix4 lRes;
	lRes << lS.x(), lS.y(), lS.z(), -lS.dot(pEye),
		lU.x(), lU.y(), lU.z(), -lU.dot(pEye),
		-lF.x(), -lF.y(), -lF.z(), lF.dot(pEye),
		0, 0, 0, 1;

	return lRes;
}

static int intersectObject(int pSI, int pSE, Matrix4f pM, Matrix4f pV, Matrix4f pP, double Px, double Py)
{
	for (int i = pSI; i < pSE; i = i + 3)
	{
		Vector4f lT1;
		lT1 << V(0, i), V(1, i), V(2, i), 1;
		Vector4f lT2;
		lT2 << V(0, i + 1), V(1, i + 1), V(2, i + 1), 1;
		Vector4f lT3;
		lT3 << V(0, i + 2), V(1, i + 2), V(2, i + 2), 1;

		Vector4f lFT1 = pP * pV * pM * lT1;
		Vector4f lFT2 = pP * pV * pM * lT2;
		Vector4f lFT3 = pP * pV * pM * lT3;

		double x1 = (double)lFT1(0);
		double y1 = (double)lFT1(1);
		double x2 = (double)lFT2(0);
		double y2 = (double)lFT2(1);
		double x3 = (double)lFT3(0);
		double y3 = (double)lFT3(1);

		if (isInsideTriangle(x1, y1, x2, y2, x3, y3, Px, Py))
		{
			//std::cout << " x1 " << x1 << " y1 " << y1 << " x2 " << x2 << " y2 " << y2 << " x3 " << x3 << " y3 " << y3 << " Px " << Px << " Py " << Py << " \n";
			return 1;
		}
	}
	return -1;
}

int getSelectedObjectNumber(double pXworld, double pYworld)
{
	int liSO = -1;

	Vector3f lEye;
	if (miCameraMode == 1)
	{
		lEye << mdCameraX + (miXPress * 0.2), mdCameraY + (miYPress * 0.2), mdCameraZ + (miZPress * 0.2);
	}
	else
	{
		double lAngleX = (30 * (miXPress % 12));     //15,24
		double lRadX = (lAngleX * PI) / 180;

		double lAngleY = (30 * (miYPress % 12));     //15,24
		double lRadY = (lAngleY * PI) / 180;

		double camX = ((miZPress * 0.5) + 2) * sin(lRadX) * cos(lRadY);
		double camY = ((miZPress * 0.5) + 2) * -sinf(lRadY);
		double camZ = ((miZPress * 0.5) + 2) * cosf(lRadX) * cosf(lRadY);

		lEye << camX, camY, camZ;
	}

	Vector3f lCenter;
	lCenter << 0, 0, 0;
	Vector3f lEyeUp;
	lEyeUp << 0, 1, 0;
	Matrix4f lCameraMatrix = getLookAtMatrix(lEye, lCenter, lEyeUp);

	Matrix4f lOrthoProjMatrix;
	lOrthoProjMatrix << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	for (int i = 1; i <= miNoOfObjects; i++)
	{
		int lSindex = getStartIndex(i);
		int lEindex = V.cols() - 1;
		if (i != miNoOfObjects)
		{
			lEindex = getStartIndex(i + 1);
		}

		//Rotation
		Matrix4f lRotationMatrix;
		double pAngleDegree = 0;
		int noOfRotation = mlistObjR.at(i - 1);
		if (noOfRotation >= 0)
		{
			for (int j = 0; j < noOfRotation;j++)
			{
				pAngleDegree = pAngleDegree + 20;
			}
		}
		else
		{
			for (int j = noOfRotation; j < 0;j++)
			{
				pAngleDegree = pAngleDegree - 20;
			}
		}
		double pAngleInRadians = (pAngleDegree * PI) / 180;
		lRotationMatrix << cos(pAngleInRadians), -sin(pAngleInRadians), 0, 0,
			sin(pAngleInRadians), cos(pAngleInRadians), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;

		//Scaling
		Matrix4f lScalingMatrix;
		int noOfScaling = mlistObjS.at(i - 1);
		double ldScale = 1;
		if (noOfScaling >= 0)
		{
			for (int j = 0; j < noOfScaling;j++)
			{
				ldScale = ldScale + (0.2 * ldScale);
			}
		}
		else
		{
			for (int j = noOfScaling; j < 0;j++)
			{
				ldScale = ldScale - (0.2 * ldScale);
			}
		}

		lScalingMatrix << ldScale, 0, 0, 0,
			0, ldScale, 0, 0,
			0, 0, ldScale, 0,
			0, 0, 0, 1;

		//Translation
		Matrix4f lTranslationMatrix;
		//std::cout << mlistObjTX.at(i) << " x " << mlistObjTY.at(i) << " -= " << i << "\n";
		lTranslationMatrix << 1, 0, 0, mlistObjTX.at(i - 1),
			0, 1, 0, mlistObjTY.at(i - 1),
			0, 0, 1, 0,
			0, 0, 0, 1;

		//Model Matrix
		Matrix4f lModelMatrix = lTranslationMatrix * lRotationMatrix * lScalingMatrix;

		if (intersectObject(lSindex, lEindex, lModelMatrix, lCameraMatrix, lOrthoProjMatrix, pXworld, pYworld) == 1)
		{
			liSO = i;
			break;
		}
	}

	return liSO;
}




void revertOldObjectColor(int pObjectNumber)
{
	if (pObjectNumber == -1)
	{
		return;
	}
	int sIndex = getStartIndex(pObjectNumber);
	if (mlistObjType.at(pObjectNumber - 1) == 1)
	{
		for (int j = 0; j < 36; j++)
		{
			C.col(sIndex + j) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
		}
	}
	else
	{
		//hardcoding as 3000 since bunny and bumpy as 1000faces*3vertices which is constant, could have use liBuTotalTriangle,liBcTotalTriangle but currently hardcoding it. 
		for (int j = 0; j < 3000; j++)
		{
			C.col(sIndex + j) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
		}
	}
	VBO_C.update(C);
}


void changeNewObjectColor(int pObjectNumber)
{
	if (pObjectNumber == -1)
	{
		return;
	}
	int sIndex = getStartIndex(pObjectNumber);
	if (mlistObjType.at(pObjectNumber - 1) == 1)
	{
		for (int j = 0; j < 36; j++)
		{
			C.col(sIndex + j) << mGreenColor(0), mGreenColor(1), mGreenColor(2);
		}
	}
	else
	{
		//hardcoding as 3000 since bunny and bumpy as 1000faces*3vertices which is constant, could have use liBuTotalTriangle,liBcTotalTriangle but currently hardcoding it. 
		for (int j = 0; j < 3000; j++)
		{
			C.col(sIndex + j) << mGreenColor(0), mGreenColor(1), mGreenColor(2);
		}
	}
	VBO_C.update(C);
}


int convertToSVGCoordinatesX(double pX)
{
	double ldSX = -1.0;
	double diff = pX - ldSX;
	return 100 * diff;
}

int convertToSVGCoordinatesY(double pY)
{
	double ldSY = 1.0;
	double diff = pY - ldSY;
	return 100 * (diff * -1);
}



void createSVG()
{
	Vector3f lEye;
	if (miCameraMode == 1)
	{
		lEye << mdCameraX + (miXPress * 0.2), mdCameraY + (miYPress * 0.2), mdCameraZ + (miZPress * 0.2);
	}
	else
	{
		double lAngleX = (30 * (miXPress % 12));     //15,24
		double lRadX = (lAngleX * PI) / 180;

		double lAngleY = (30 * (miYPress % 12));     //15,24
		double lRadY = (lAngleY * PI) / 180;

		double camX = ((miZPress * 0.5) + 2) * sin(lRadX) * cos(lRadY);
		double camY = ((miZPress * 0.5) + 2) * -sinf(lRadY);
		double camZ = ((miZPress * 0.5) + 2) * cosf(lRadX) * cosf(lRadY);

		lEye << camX, camY, camZ;
	}

	//	lEye << 0, 0, 1;
	Vector3f lCenter;
	lCenter << 0, 0, 0;
	Vector3f lEyeUp;
	lEyeUp << 0, 1, 0;
	Matrix4f cameraMatrix = getLookAtMatrix(lEye, lCenter, lEyeUp);

	//	std::cout << width << " " << height;
	if (miHeight == 0)
	{
		miHeight = 1;
	}
	double ldAspectrRatio = (double)miWidth / miHeight;
	double ldNear = 0.1;
	double ldFar = 100.0;
	double ldFOVdegrees = 60.0;

	Matrix4f lProjectionMatrix = getPerspectiveMatrix<Matrix4f::Scalar>(
		ldFOVdegrees,
		ldAspectrRatio,
		ldNear,
		ldFar
		);

	Matrix4f lOrthoProjMatrix;
	lOrthoProjMatrix << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;



	//---------------------------------------------------
	string lStart = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?><svg width = \"200\" height = \"200\" xmlns = \"http://www.w3.org/2000/svg\" ><rect x = \"0\" y = \"0\" width = \"200\" height = \"200\" fill = \"rgb(76, 76, 76)\" />";
	string lTriangleTemplateStart = "<polygon points=\"";
	string lTTM = "\"";
	//string lTriangleTemplateEnd = " stroke=\"rgb(255, 51, 153)\" stroke-width=\"2\" fill=\"rgb(255, 51, 153)\" />";

	string lTriangleTemplateEnd1 = " stroke=\"rgb(";
	string lTriangleTemplateEnd3 = ")\" stroke-width=\"2\" fill=\"rgb(";
	
	string lTriangleTemplateEnd2 = ")\" />";




	string lend = "</svg>";
	string lfinalstring = lStart;

	//----
	for (int i = 0; i < miNoOfObjects; i++)
	{
		if (mlistObjDeleted.at(i) == 1)
		{
			continue;
		}
		int lSindex = getStartIndex(i + 1);
		int lEindex;
		if (mlistObjType.at(i) == 1)
		{
			lEindex = lSindex + 36;
		}
		else if (mlistObjDispType.at(i) == 2)
		{
			lEindex = lSindex + 3000;
		}
		else
		{
			lEindex = lSindex + 3000;
		}

		Matrix4f lRotationMatrix;
		double pAngleDegree = 0;
		int noOfRotation = mlistObjR.at(i);
		if (noOfRotation >= 0)
		{
			for (int j = 0; j < noOfRotation;j++)
			{
				pAngleDegree = pAngleDegree + 20;
			}
		}
		else
		{
			for (int j = noOfRotation; j < 0;j++)
			{
				pAngleDegree = pAngleDegree - 20;
			}
		}
		double pAngleInRadians = (pAngleDegree * PI) / 180;
		lRotationMatrix << cos(pAngleInRadians), -sin(pAngleInRadians), 0, 0,
			sin(pAngleInRadians), cos(pAngleInRadians), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;

		//Scaling
		Matrix4f lScalingMatrix;
		int noOfScaling = mlistObjS.at(i);
		double ldScale = 1;
		if (noOfScaling >= 0)
		{
			for (int j = 0; j < noOfScaling;j++)
			{
				ldScale = ldScale + (0.2 * ldScale);
			}
		}
		else
		{
			for (int j = noOfScaling; j < 0;j++)
			{
				ldScale = ldScale - (0.2 * ldScale);
			}
		}

		lScalingMatrix << ldScale, 0, 0, 0,
			0, ldScale, 0, 0,
			0, 0, ldScale, 0,
			0, 0, 0, 1;

		//Translation
		Matrix4f lTranslationMatrix;
		//	std::cout << mlistObjTX.at(i) << " x " << mlistObjTY.at(i) << " -= " << i << "\n";

		lTranslationMatrix << 1, 0, 0, mlistObjTX.at(i),
			0, 1, 0, mlistObjTY.at(i),
			0, 0, 1, 0,
			0, 0, 0, 1;

		//Model Matrix
		Matrix4f lModelMatrix = lTranslationMatrix * lRotationMatrix * lScalingMatrix; // * Translation


		for (int k = lSindex; k < lEindex; k = k + 3)
		{
			Vector4f lT1;
			lT1 << V.col(k)(0), V.col(k)(1), V.col(k)(2), 1;

			//Vector4f lFT1 = lT1;
			Vector4f lFT1 = lOrthoProjMatrix * cameraMatrix * lModelMatrix *  lT1;

			Vector4f lT2;
			lT2 << V.col(k + 1)(0), V.col(k + 1)(1), V.col(k + 1)(2), 1;

			//Vector4f lFT2 = lT2;
			Vector4f lFT2 = lOrthoProjMatrix * cameraMatrix * lModelMatrix * lT2;

			Vector4f lT3;
			lT3 << V.col(k + 2)(0), V.col(k + 2)(1), V.col(k + 2)(2), 1;

			//Vector4f lFT3 =  lT3;
			Vector4f lFT3 = lOrthoProjMatrix * cameraMatrix * lModelMatrix  * lT3;


			string lcValues = "";

			Vector3f lightColor;
			lightColor << 1.0, 1.0, 1.0;
			Vector3f lightPos;
			lightPos << 0.4, 0.3, 0.5;
			Vector4f fragPos1 = lModelMatrix * lFT3;

			Vector3f fragPos;
			fragPos << fragPos1(0), fragPos1(1), fragPos1(2);
			Vector3f ambient;
			ambient << 0.2, 0.2, 0.2; //vec3 ambient = ambientStrength * lightColor;
			Vector4f Normal;
			Normal << NF.col(k)(0), NF.col(k)(0), NF.col(k)(0), 1;
			Vector4f bNormal = lModelMatrix.inverse().transpose() * Normal;
			Vector3f aNormal;
			aNormal << bNormal(0), bNormal(0), bNormal(0);
			Vector3f lNorm = aNormal.normalized();
			Vector3f lightDir = (lightPos - fragPos).normalized();
			float diff = max(lNorm.dot(lightDir), 0.0f);
			Vector3f diffuse;
			diffuse << diff * lightColor;
			Vector3f viewDir = (lEye - fragPos).normalized();
			Vector3f reflectDir = (-1 * lightDir) - (2.0 * lNorm.dot(lightDir)) * lNorm;
			float spec = pow(max(viewDir.dot(reflectDir), 0.0f), 32);
			Vector3f specular;
			specular << 0.5 * spec * lightColor;
			Vector3f result1 = (ambient + diffuse + specular);
			Vector3f result;
			result << result1(0) * mLPinkColor(0), result1(1)* mLPinkColor(1), result1(2)* mLPinkColor(2);
			result = result * 255;
			//std::cout << result << "\n";

			lcValues = lcValues + to_string((int)result(0)) + "," + to_string((int)result(1)) + "," + to_string((int)result(2));


			string lValues = "";

			int x = convertToSVGCoordinatesX(lFT1(0));
			int y = convertToSVGCoordinatesY(lFT1(1));
			lValues = lValues + to_string(x) + "," + to_string(y) + ",";

			x = convertToSVGCoordinatesX(lFT2(0));
			y = convertToSVGCoordinatesY(lFT2(1));
			lValues = lValues + to_string(x) + "," + to_string(y) + ",";

			x = convertToSVGCoordinatesX(lFT3(0));
			y = convertToSVGCoordinatesY(lFT3(1));
			lValues = lValues + to_string(x) + "," + to_string(y);

			string lfinaltriangle = lTriangleTemplateStart + lValues + lTTM + lTriangleTemplateEnd1 + lcValues  + lTriangleTemplateEnd3 + lcValues + lTriangleTemplateEnd2;
			lfinalstring = lfinalstring + lfinaltriangle;
		}	
	}


	lfinalstring = lfinalstring + lend;
	//std::cout << lfinalstring << "\n";
	ofstream file("output.svg");
	file << lfinalstring;
	file.close();

}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates
	double xcworld = ((xpos / double(width)) * 2) - 1;
	double ycworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw

	if (miCurrentMode == 2 && mbIsPressed == 1 && miSelectedObject != -1)
	{
		double mdChangeTriangleX = xcworld - mdSelectedObjectX;
		double mdChangeTriangleY = ycworld - mdSelectedObjectY;

		//	std::cout << mdChangeTriangleX << " " << mdChangeTriangleY << "\n";
		tempno = 0;
		tempno = (double)mlistObjTX.at(miSelectedObject - 1);
		tempno = tempno + mdChangeTriangleX;
		mlistObjTX.at(miSelectedObject - 1) = (float)tempno;

		tempno = 0;
		tempno = (double)mlistObjTY.at(miSelectedObject - 1);
		tempno = tempno + mdChangeTriangleY;
		mlistObjTY.at(miSelectedObject - 1) = (float)tempno;

		mdSelectedObjectX = xcworld;
		mdSelectedObjectY = ycworld;
	}

}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the position of the mouse in the window
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	//std::cout << xpos << " " << ypos << "\n";
	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	//std::cout << width << " " << height << "\n";
	// Convert screen position to world coordinates
	double xworld = ((xpos / double(width)) * 2) - 1;
	double yworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw
	//std::cout << xworld << " " << yworld << "\n";

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		mbIsPressed = 1;

		if (miCurrentMode == 2)
		{
			int liOldObject = miSelectedObject;
			miSelectedObject = getSelectedObjectNumber(xworld, yworld);
			revertOldObjectColor(liOldObject);
			changeNewObjectColor(miSelectedObject);
			if (miSelectedObject != -1)
			{
				mdSelectedObjectX = xworld;
				mdSelectedObjectY = yworld;
			}
			else
			{
				mdSelectedObjectX = 0;
				mdSelectedObjectY = 0;
			}
		}
	}

	static int lsiOldState = GLFW_RELEASE;
	int liNewState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if (liNewState == GLFW_RELEASE && lsiOldState == GLFW_PRESS) {
		mbIsPressed = 0;
	}
	lsiOldState = liNewState;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_I:
			std::cout << "I is pressed" << "\n";
			miCurrentMode = 1;
			break;
		case GLFW_KEY_1:
			std::cout << "1 is pressed" << "\n";
			if (miCurrentMode == 1)
			{
				addCube();
			}
			break;
		case GLFW_KEY_2:
			std::cout << "2 is pressed" << "\n";
			if (miCurrentMode == 1)
			{
				addBunny();
			}
			break;
		case GLFW_KEY_3:
			std::cout << "3 is pressed" << "\n";
			if (miCurrentMode == 1)
			{
				addBumpyCube();
			}
			break;
		case GLFW_KEY_O:
			std::cout << "O is pressed" << "\n";
			miCurrentMode = 2;
			break;
		case GLFW_KEY_H:
			std::cout << "H is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			tempno = mlistObjR.at(miSelectedObject - 1);
			tempno++;
			mlistObjR.at(miSelectedObject - 1) = tempno;
			break;
		case GLFW_KEY_J:
			std::cout << "J is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			tempno = mlistObjR.at(miSelectedObject - 1);
			tempno--;
			mlistObjR.at(miSelectedObject - 1) = tempno;
			//	for (auto v : mlistObjR)
			//		std::cout << v << " ";
			break;
		case GLFW_KEY_K:
			std::cout << "K is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			tempno = mlistObjS.at(miSelectedObject - 1);
			tempno++;
			mlistObjS.at(miSelectedObject - 1) = tempno;
			break;
		case GLFW_KEY_L:
			std::cout << "L is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			tempno = mlistObjS.at(miSelectedObject - 1);
			tempno--;
			mlistObjS.at(miSelectedObject - 1) = tempno;
			break;
		case GLFW_KEY_P:
			std::cout << "P is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			mlistObjDeleted.at(miSelectedObject - 1) = 1;
			break;
		case GLFW_KEY_C:
			std::cout << "C is pressed" << "\n";
			miCameraMode = 2;
			break;
		case GLFW_KEY_Z:
			std::cout << "Z is pressed" << "\n";
			miProjectionMode = 1;
			break;
		case GLFW_KEY_X:
			std::cout << "X is pressed" << "\n";
			miProjectionMode = 2;
			break;
		case GLFW_KEY_V:
			miXPress = 0;
			miYPress = 0;
			miZPress = 0;
			miCameraMode = 1;
			break;
		case GLFW_KEY_W:
			std::cout << "W is pressed" << "\n";
			miYPress = miYPress + 1;
			break;
		case GLFW_KEY_S:
			std::cout << "S is pressed" << "\n";
			miYPress = miYPress - 1;
			break;
		case GLFW_KEY_A:
			std::cout << "A is pressed" << "\n";
			miXPress = miXPress - 1;
			break;
		case GLFW_KEY_D:
			std::cout << "D is pressed" << "\n";
			miXPress = miXPress + 1;
			break;
		case GLFW_KEY_Q:
			std::cout << "Q is pressed" << "\n";
			miZPress = miZPress + 1;
			break;
		case GLFW_KEY_E:
			std::cout << "E is pressed" << "\n";
			miZPress = miZPress - 1;
			break;
		case GLFW_KEY_R:
			std::cout << "R is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			mlistObjDispType.at(miSelectedObject - 1) = 1;
			break;
		case GLFW_KEY_T:
			std::cout << "T is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			mlistObjDispType.at(miSelectedObject - 1) = 2;
			break;
		case GLFW_KEY_Y:
			std::cout << "Y is pressed" << "\n";
			if (miCurrentMode != 2 || miSelectedObject == -1)
			{
				return;
			}
			mlistObjDispType.at(miSelectedObject - 1) = 3;
			break;
		case GLFW_KEY_G:
			std::cout << "G is pressed" << "\n";
			createSVG();
			break;
		case GLFW_KEY_B:
			std::cout << "B is pressed" << "\n";
			miCurrentMode = 4;
			break;
		default:
			break;
		}
	}
}


double bezier(double pA, double pB, double pC, double pD, double pT)
{
	//std::cout << pA << pB << pC << pD << pT << "\n";
	double lS = 1 - pT;
	double lAB = pA * lS + pB * pT;
	double lBC = pB * lS + pC * pT;
	double lCD = pC * lS + pD * pT;
	double lABC = lAB * lS + lCD * pT;
	double lBCD = lBC * lS + lCD * pT;
	return lABC * lS + lBCD * pT;
}

double getBCXYValues(int XY, double pI)
{
	if (XY == 1)
	{
		double x = bezier(0, 0, 1, 1, pI);
		return x;
	}
	else if (XY == 2)
	{
		double y = bezier(0, -1, -1, 0, pI);
		return y;
	}
}

int main(void)
{
	GLFWwindow* window;

	// Initialize the library
	if (!glfwInit())
		return -1;

	// Activate supersampling
	glfwWindowHint(GLFW_SAMPLES, 8);

	// Ensure that we get at least a 3.2 context
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

	// On apple we have to load a core profile with forward compatibility
#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

	// Create a windowed mode window and its OpenGL context
	window = glfwCreateWindow(640, 480, "Rupesh World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	// Make the window's context current
	glfwMakeContextCurrent(window);

#ifndef __APPLE__
	glewExperimental = true;
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}
	glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
	fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
#endif

	int major, minor, rev;
	major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
	minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
	rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
	printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
	printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
	printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));


	// Initialize the VAO
	// A Vertex Array Object (or VAO) is an object that describes how the vertex
	// attributes are stored in a Vertex Buffer Object (or VBO). This means that
	// the VAO is not the actual object storing the vertex data,
	// but the descriptor of the vertex data.
	VertexArrayObject VAO;
	VAO.init();
	VAO.bind();

	// Initialize the VBO with the vertices data
	// A VBO is a data container that lives in the GPU memory
	VBO.init();

	// Second VBO for colors
	VBO_C.init();
	VBO_NF.init();
	VBO_NP.init();

	VBO.update(V);
	VBO_C.update(C);
	VBO_NF.update(NF);
	VBO_NP.update(NP);

	// Initialize the OpenGL Program
	// A program controls the OpenGL pipeline and it must contains
	// at least a vertex shader and a fragment shader to be valid
	Program program;
	const GLchar* vertex_shader =
		"#version 150 core\n"
		"in vec3 position;"
		"in vec3 color;"
		"in vec3 normal1;"
		"in vec3 normal2;"
		"uniform mat4 view;"
		"uniform mat4 projection;"
		"uniform mat4 model;"
		"uniform float shading;"
		"uniform float applyShading;"
		"uniform vec3 lightPos;"
		"uniform vec3 viewPos;"
		"uniform vec3 lightColor;"
		"out vec3 f_color;"
		"vec3 normal;"
		"vec3 fragPos;"
		"vec3 aNormal;"
		"float ambientStrength = 0.2;"
		"float specularStrength = 0.5;"
		"void main()"
		"{"
		"if (shading == 2.0)"
		"{"
		"	normal = normal2;"
		"}"
		"else if (shading == 1.0)"
		"{"
		"	normal = normal1;"
		"}"
		"    gl_Position =  projection * view * model * vec4(position, 1.0);"
		"	 fragPos = vec3( model * vec4(position, 1.0) );"
		"	 vec3 ambient = ambientStrength * lightColor;"
		"	 aNormal = mat3(transpose(inverse(model))) * normal;"
		"	 vec3 norm = normalize(aNormal);"
		"	 vec3 lightDir = normalize(lightPos - fragPos);"
		"	 float diff = max(dot(norm, lightDir), 0.0);"
		"	 vec3 diffuse = diff * lightColor;"
		"    vec3 viewDir = normalize(viewPos - fragPos);"
		"	 vec3 reflectDir = reflect(-lightDir, norm);"
		"	 float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);"
		"	 vec3 specular = specularStrength * spec * lightColor;"
		"	 vec3 result = color;"
		"	if (applyShading != 0.0)"
		"	{"
		"		result = (ambient + diffuse + specular) * color;"
		"	}"
		"    f_color = result;"
		"}";
	const GLchar* fragment_shader =
		"#version 150 core\n"
		"in vec3 f_color;"
		"out vec4 outColor;"
		"void main()"
		"{"
		"    outColor = vec4(f_color, 1.0);"
		"}";

	// Compile the two shaders and upload the binary to the GPU
	// Note that we have to explicitly specify that the output "slot" called outColor
	// is the one that we want in the fragment buffer (and thus on screen)
	program.init(vertex_shader, fragment_shader, "outColor");
	program.bind();

	// The vertex shader wants the position of the vertices as an input.
	// The following line connects the VBO we defined above with the position "slot"
	// in the vertex shader
	program.bindVertexAttribArray("position", VBO);
	program.bindVertexAttribArray("color", VBO_C);
	program.bindVertexAttribArray("normal1", VBO_NF);
	program.bindVertexAttribArray("normal2", VBO_NP);

	// Register the keyboard callback
	glfwSetKeyCallback(window, key_callback);

	// Register the mouse callback
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);

	// Update viewport
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	Vector3f  lLightPos;
	lLightPos << 0.4, 0.3, 0.5;
	glUniform3fv(program.uniform("lightPos"), 1, lLightPos.data());

	Vector3f  lLightColor;
	lLightColor << 1.0, 1.0, 1.0;
	glUniform3fv(program.uniform("lightColor"), 1, lLightColor.data());



	// Save the current time --- it will be used to dynamically change the triangle color
	auto t_start = std::chrono::high_resolution_clock::now();
	int timeCounter = 0;
	float oldtime = 0;

		// Loop until the user closes the window
	while (!glfwWindowShouldClose(window))
	{
		// Bind your VAO (not necessary if you have only one)
		VAO.bind();

		// Bind your program
		program.bind();

		// Clear the framebuffer
		glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Enable depth test
		glEnable(GL_DEPTH_TEST);

		Vector3f lEye;
		if (miCameraMode == 1)
		{
			lEye << mdCameraX + (miXPress * 0.2), mdCameraY + (miYPress * 0.2), mdCameraZ + (miZPress * 0.2);
		}
		else
		{
			double lAngleX = (30 * (miXPress % 12));     //15,24
			double lRadX = (lAngleX * PI) / 180;

			double lAngleY = (30 * (miYPress % 12));     //15,24
			double lRadY = (lAngleY * PI) / 180;

			double camX = ((miZPress * 0.5) + 2) * sin(lRadX) * cos(lRadY);
			double camY = ((miZPress * 0.5) + 2) * -sinf(lRadY);
			double camZ = ((miZPress * 0.5) + 2) * cosf(lRadX) * cosf(lRadY);

			lEye << camX, camY, camZ;
		}

		//	lEye << 0, 0, 1;
		Vector3f lCenter;
		lCenter << 0, 0, 0;
		Vector3f lEyeUp;
		lEyeUp << 0, 1, 0;
		Matrix4f cameraMatrix = getLookAtMatrix(lEye, lCenter, lEyeUp);

		glUniformMatrix4fv(program.uniform("view"), 1, GL_FALSE, cameraMatrix.data());

		glUniform3fv(program.uniform("viewPos"), 1, lEye.data());


		int liWidth, liHeight;
		glfwGetWindowSize(window, &liWidth, &liHeight);
		//	std::cout << width << " " << height;
		if (liHeight == 0)
		{
			liHeight = 1;
		}
		double ldAspectrRatio = (double)liWidth / liHeight;
		double ldNear = 0.1;
		double ldFar = 100.0;
		double ldFOVdegrees = 60.0;

		Matrix4f lProjectionMatrix = getPerspectiveMatrix<Matrix4f::Scalar>(
			ldFOVdegrees,
			ldAspectrRatio,
			ldNear,
			ldFar
			);

		Matrix4f lOrthoProjMatrix;
		lOrthoProjMatrix << 1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1;


		if (miProjectionMode == 1)
		{
			glUniformMatrix4fv(program.uniform("projection"), 1, GL_FALSE, lOrthoProjMatrix.data());
		}
		else
		{
			glUniformMatrix4fv(program.uniform("projection"), 1, GL_FALSE, lProjectionMatrix.data());
		}

		//std::cout << miNoOfObjects << "\n";
		for (int i = 0; i < miNoOfObjects; i++)
		{
			if (mlistObjDeleted.at(i) == 1)
			{
				continue;
			}
			int lSindex = getStartIndex(i + 1);
			int lEindex;
			if (mlistObjType.at(i) == 1)
			{
				lEindex = lSindex + 36;
			}
			else if (mlistObjDispType.at(i) == 2)
			{
				lEindex = lSindex + 3000;
			}
			else
			{
				lEindex = lSindex + 3000;
			}

			//Rotation
			Matrix4f lRotationMatrix;
			double pAngleDegree = 0;
			int noOfRotation = mlistObjR.at(i);
			if (noOfRotation >= 0)
			{
				for (int j = 0; j < noOfRotation;j++)
				{
					pAngleDegree = pAngleDegree + 20;
				}
			}
			else
			{
				for (int j = noOfRotation; j < 0;j++)
				{
					pAngleDegree = pAngleDegree - 20;
				}
			}
			double pAngleInRadians = (pAngleDegree * PI) / 180;
			lRotationMatrix << cos(pAngleInRadians), -sin(pAngleInRadians), 0, 0,
				sin(pAngleInRadians), cos(pAngleInRadians), 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1;

			//Scaling
			Matrix4f lScalingMatrix;
			int noOfScaling = mlistObjS.at(i);
			double ldScale = 1;
			if (noOfScaling >= 0)
			{
				for (int j = 0; j < noOfScaling;j++)
				{
					ldScale = ldScale + (0.2 * ldScale);
				}
			}
			else
			{
				for (int j = noOfScaling; j < 0;j++)
				{
					ldScale = ldScale - (0.2 * ldScale);
				}
			}

			lScalingMatrix << ldScale, 0, 0, 0,
				0, ldScale, 0, 0,
				0, 0, ldScale, 0,
				0, 0, 0, 1;

			// Set the uniform value depending on the time difference
			auto t_now = std::chrono::high_resolution_clock::now();
			float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
		//std:cout << " x " << getBCXYValues(1, timeCounter * 0.2) << " y " << getBCXYValues(2, timeCounter * 0.2) << " i " << (timeCounter * 0.2) << " \n";


			//Translation
			Matrix4f lTranslationMatrix;
			//	std::cout << mlistObjTX.at(i) << " x " << mlistObjTY.at(i) << " -= " << i << "\n";

			if (miCurrentMode == 4)
			{
				if (time - oldtime > 2)
				{
					timeCounter++;
					oldtime = time;
				}
				if (timeCounter > 5)
				{
					timeCounter = 0;
				}
				lTranslationMatrix << 1, 0, 0, mlistObjTX.at(i) + getBCXYValues(1, (timeCounter * 0.2)),
					0, 1, 0, mlistObjTY.at(i) + getBCXYValues(2, (timeCounter * 0.2)),
					0, 0, 1, 0,
					0, 0, 0, 1;
			}
			else
			{
				lTranslationMatrix << 1, 0, 0, mlistObjTX.at(i),
					0, 1, 0, mlistObjTY.at(i),
					0, 0, 1, 0,
					0, 0, 0, 1;
			}

			

			//Model Matrix
			Matrix4f lModelMatrix = lTranslationMatrix * lRotationMatrix * lScalingMatrix; // * Translation
			glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, lModelMatrix.data());

			if (mlistObjDispType.at(i) == 1)
			{
				glUniform1f(program.uniform("shading"), 1.0);
				glUniform1f(program.uniform("applyShading"), 0.0);

				glLineWidth(3.0);
				glDrawArrays(GL_LINE_STRIP, lSindex, lEindex);
				
				
			}
			if (mlistObjDispType.at(i) == 2)
			{
				glUniform1f(program.uniform("shading"), 1.0);
				glUniform1f(program.uniform("applyShading"), 1.0);
				glDrawArrays(GL_TRIANGLES, lSindex, lEindex);
				//				glLineWidth(3.0);
				//				glDrawArrays(GL_LINE_STRIP, lSindex, lEindex);
			}
			else if (mlistObjDispType.at(i) == 3)
			{
				glUniform1f(program.uniform("shading"), 2.0);
				glUniform1f(program.uniform("applyShading"), 1.0);
			//	std::cout << i << " " << lEindex << "\n";
				glDrawArrays(GL_TRIANGLES, lSindex, lEindex);
			}
		}

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();
	}

	// Deallocate opengl memory
	program.free();
	VAO.free();
	VBO.free();

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}
