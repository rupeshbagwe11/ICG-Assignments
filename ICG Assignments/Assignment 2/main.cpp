// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include<cmath>
#include <limits>
#include <fstream>


#include <Eigen/Core>
#include <Eigen/Geometry>


using namespace std;
using namespace Eigen;

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif

// Timer
#include <chrono>

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_C;

// Contains the vertex positions
Eigen::MatrixXf V(2, 1);
Eigen::MatrixXf V_bkp(2, 1);

// Contains the per-vertex color
Eigen::MatrixXf C(3, 1);
Eigen::Matrix3f OC(3, 3);

int miStepCounter = 0;
// i = 1, o = 2, p = 3, c = 4, r = 5, b = 6
int miCurrentMode = 1;
int miSelectedTriangle = -1;
int miClosestVertex = -1;
int mbIsPressed = 0;

int miRotationMode = -1;
int miScalingMode = -1;

double mdSelectedTriangleX = -2.0;
double mdSelectedTriangleY = -2.0;

double mdChangeTriangleX = 0.0;
double mdChangeTriangleY = 0.0;

double mdScaleTriangleX = 0.0;
double mdScaleTriangleY = 0.0;

double mdLR = 0;
double mdUD = 0;
double mdZoom = 0;

double pi = 3.14159;

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


void changeVertexColor(int pVertexNo, Vector3d pColor)
{
	if (pVertexNo == -1)
	{
		return;
	}
	C.col(pVertexNo) << pColor(0), pColor(1), pColor(2);
}

void changeTriangleColorMul(int pTriangleNo, Matrix3f pColor)
{
	int liTriangleColNo = (pTriangleNo * 3) - 1;
	int j = 0;
	for (int i = 2; i >= 0; i--)
	{
		int liIndex = liTriangleColNo - i;
		Eigen::MatrixXd TC(3, 1);
		TC << pColor.col(j)(0), pColor.col(j)(1), pColor.col(j)(2);
		j++;
		changeVertexColor(liIndex, TC);
	}
}

void changeTriangleColor(int pTriangleNo, Vector3d pColor)
{
	int liTriangleColNo = (pTriangleNo * 3) - 1;
	for (int i = 2; i >= 0; i--)
	{
		int liIndex = liTriangleColNo - i;
		changeVertexColor(liIndex, pColor);
	}
}

void saveOldTriangleColor(int pTriangleNo)
{
	int liTriangleColNo = (pTriangleNo * 3) - 1;
	int liIndex = liTriangleColNo - 2;
	//std::cout << C << "\n";
	//std::cout << liIndex << "\n";

	OC.col(0)(0) = C.col(liIndex)(0);
	OC.col(0)(1) = C.col(liIndex)(1);
	OC.col(0)(2) = C.col(liIndex)(2);
	OC.col(1)(0) = C.col(liIndex+1)(0);
	OC.col(1)(1) = C.col(liIndex+1)(1);
	OC.col(1)(2) = C.col(liIndex+1)(2);
	OC.col(2)(0) = C.col(liIndex+2)(0);
	OC.col(2)(1) = C.col(liIndex+2)(1);
	OC.col(2)(2) = C.col(liIndex+2)(2);

	//std::cout << C.col(liIndex)(0) << "\n";
	//std::cout << liIndex << "\n";
	//std::cout << OC << "\n";
}

void selectTriangleColor(int pOldTriangleNo, int pNewTriangleNo)
{

	if (pOldTriangleNo != pNewTriangleNo)
	{
		//std::cout << pOldTriangleNo << "  " << pNewTriangleNo << "\n";
		if (pOldTriangleNo != -1)
		{
			changeTriangleColorMul(pOldTriangleNo, OC);
			//changeTriangleColor1(pOldTriangleNo, mLPinkColor);
		}

		if (pNewTriangleNo != -1)
		{
			saveOldTriangleColor(pNewTriangleNo);
			changeTriangleColor(pNewTriangleNo, mLBlueColor);
		}
	}
}

void refreshMembers()
{
	if (miCurrentMode != 2 && miSelectedTriangle != -1)
	{
		int liOldSelectedTriangle = miSelectedTriangle;
		miSelectedTriangle = -1;
		selectTriangleColor(liOldSelectedTriangle, miSelectedTriangle);
	}
	mdSelectedTriangleX = -2.0;
	mdSelectedTriangleY = -2.0;
	mdChangeTriangleX = 0.0;
	mdChangeTriangleY = 0.0;
	mdScaleTriangleX = 0.0;
	mdScaleTriangleY = 0.0;
	miRotationMode = -1;
	miScalingMode = -1;
	if (miCurrentMode != 4 && miClosestVertex != -1)
	{
		miClosestVertex = -1;
	}
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void removeColumn(Eigen::MatrixXf& pMatrix, unsigned int pColToRemove)
{
	unsigned int liNumRows = pMatrix.rows();
	unsigned int liNumCols = pMatrix.cols() - 1;

	if (pColToRemove < liNumCols)
		pMatrix.block(0, pColToRemove, liNumRows, liNumCols - pColToRemove) = pMatrix.block(0, pColToRemove + 1, liNumRows, liNumCols - pColToRemove);

	pMatrix.conservativeResize(liNumRows, liNumCols);
}

static void deleteTriangle(int pTriangleNo)
{
	refreshMembers();
	//std::cout << V << "\n";
	removeColumn(V, (pTriangleNo * 3) - 1);
	removeColumn(C, (pTriangleNo * 3) - 1);
	miStepCounter--;
	removeColumn(V, (pTriangleNo * 3) - 2);
	removeColumn(C, (pTriangleNo * 3) - 2);
	miStepCounter--;
	removeColumn(V, (pTriangleNo * 3) - 3);
	removeColumn(C, (pTriangleNo * 3) - 3);
	miStepCounter--;
	//std::cout << "-----------------" << pTriangleNo << "--------------------- " << "\n";
	//std::cout << V << "\n";
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

	/* A1+A2+A3 = A with 0.002 calculateError */
	return (ldAT <= ldAST1 + ldAST2 + ldAST3 + 0.002) && (ldAT >= ldAST1 + ldAST2 + ldAST3 - 0.002);
}

static double getDistance(double pOx, double pOy, double pDx, double pDy)
{
	return sqrt(((pDx - pOx) * (pDx - pOx)) + ((pDy - pOy) * (pDy - pOy)));
}

static int getClosestVertex(double Px, double Py)
{
	double ldTempClosestValue = 10000000000000.0;
	int liVertexNo = -1;

	for (int i = 0; i < miStepCounter; i++)
	{
		double liDistance = getDistance(V.col(i)(0), V.col(i)(1), Px, Py);

		if (liDistance < ldTempClosestValue)
		{
			ldTempClosestValue = liDistance;
			liVertexNo = i;
		}
	}
	return liVertexNo;
}


static int getTriangleNumber(double Px, double Py)
{
	int liNoOfTriangles = V.cols() / 3;

	for (int i = 0; i < liNoOfTriangles; i++)
	{
		int liCurrTriangleNumber = liNoOfTriangles - i;
		int liCurrTriangleCoord = (liCurrTriangleNumber * 3) - 1;
		double x1 = V(0, liCurrTriangleCoord);
		double y1 = V(1, liCurrTriangleCoord);
		double x2 = V(0, liCurrTriangleCoord - 1);
		double y2 = V(1, liCurrTriangleCoord - 1);
		double x3 = V(0, liCurrTriangleCoord - 2);
		double y3 = V(1, liCurrTriangleCoord - 2);

		if (isInsideTriangle(x1, y1, x2, y2, x3, y3, Px, Py))
		{
			return liCurrTriangleNumber;
		}
	}

	return -1;
}

void updateTranslationTriangle(int pTriangleNo)
{
	int liTriangleColNo = (pTriangleNo * 3) - 1;
	Matrix4f lTranslationMatrix;
	lTranslationMatrix << 1, 0, 0, mdChangeTriangleX,
		0, 1, 0, mdChangeTriangleY,
		0, 0, 1, 0,
		0, 0, 0, 1;

	for (int i = 2; i >= 0; i--)
	{
		Vector4f lTemp4VVector;
		int liIndex = liTriangleColNo - i;
		lTemp4VVector << V.col(liIndex)(0), V.col(liIndex)(1), 0, 1;
		Vector4f lFinalVectorValue = lTranslationMatrix * lTemp4VVector;
		V.col(liIndex) << lFinalVectorValue(0), lFinalVectorValue(1);
	}

}


void printCentroid(int pNo)
{
	int liTriangleColNo = (miSelectedTriangle * 3) - 1;
	double lX = (V.col(liTriangleColNo - 2)(0) + V.col(liTriangleColNo - 1)(0) + V.col(liTriangleColNo)(0)) / 3.0;
	double lY = (V.col(liTriangleColNo - 2)(1) + V.col(liTriangleColNo - 1)(1) + V.col(liTriangleColNo)(1)) / 3.0;
	std::cout << std::fixed;
	std::cout << pNo << " :  "  << lX << "  " << lY << "\n";
	//std::cout << pNo << " :  " << std::cout.precision(5) << lX << "  " << std::cout.precision(5) << lY << "\n";

}


void updateScalingTriangle()
{
	if (miSelectedTriangle == -1)
	{
		refreshMembers();
		return;
	}
	//	printCentroid();
	int liTriangleColNo = (miSelectedTriangle * 3) - 1;
	double lCentroidX = (V.col(liTriangleColNo - 2)(0) + V.col(liTriangleColNo - 1)(0) + V.col(liTriangleColNo)(0)) / 3.0;
	double lCentroidY = (V.col(liTriangleColNo - 2)(1) + V.col(liTriangleColNo - 1)(1) + V.col(liTriangleColNo)(1)) / 3.0;

//	printCentroid(1);

	mdChangeTriangleX = 0 - lCentroidX;
	mdChangeTriangleY = 0 - lCentroidY;

	updateTranslationTriangle(miSelectedTriangle);

//	printCentroid(2);

	mdScaleTriangleX = 1.25;
	mdScaleTriangleY = 1.25;

	if (miScalingMode == 2)
	{
		mdScaleTriangleX = 0.75;
		mdScaleTriangleY = 0.75;
	}

	Matrix4f lScalingMatrix;
	lScalingMatrix << mdScaleTriangleX, 0, 0, 0,
		0, mdScaleTriangleY, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	for (int i = 2; i >= 0; i--)
	{
		Vector4f lTemp4VVector;
		int liIndex = liTriangleColNo - i;
		lTemp4VVector << V.col(liIndex)(0), V.col(liIndex)(1), 0, 1;
		Vector4f lFinalVectorValue = lScalingMatrix * lTemp4VVector;
		V.col(liIndex) << lFinalVectorValue(0), lFinalVectorValue(1);
	}

	//printCentroid(3);

	//double lNewCentroidX = (V.col(liTriangleColNo - 2)(0) + V.col(liTriangleColNo - 1)(0) + V.col(liTriangleColNo)(0)) / 3.0;
	//double lNewCentroidY = (V.col(liTriangleColNo - 2)(1) + V.col(liTriangleColNo - 1)(1) + V.col(liTriangleColNo)(1)) / 3.0;

	mdChangeTriangleX = lCentroidX - 0;
	mdChangeTriangleY = lCentroidY - 0;

	updateTranslationTriangle(miSelectedTriangle);

	//printCentroid(4);

	refreshMembers();
}

void updateRotationTriangle()
{
	if (miSelectedTriangle == -1)
	{
		refreshMembers();
		return;
	}

	int liTriangleColNo = (miSelectedTriangle * 3) - 1;
	double lCentroidX = (V.col(liTriangleColNo - 2)(0) + V.col(liTriangleColNo - 1)(0) + V.col(liTriangleColNo)(0)) / 3.0;
	double lCentroidY = (V.col(liTriangleColNo - 2)(1) + V.col(liTriangleColNo - 1)(1) + V.col(liTriangleColNo)(1)) / 3.0;

	//printCentroid(0);

	mdChangeTriangleX = 0 - lCentroidX;
	mdChangeTriangleY = 0 - lCentroidY;

	updateTranslationTriangle(miSelectedTriangle);

	//printCentroid(1);

	double pAngleInRadians = (10 * pi) / 180;
	if (miRotationMode == 2)
	{
		pAngleInRadians = (-10 * pi) / 180;
	}
	Matrix4f lRotationMatrix;
	lRotationMatrix << cos(pAngleInRadians), -sin(pAngleInRadians), 0, 0,
		sin(pAngleInRadians), cos(pAngleInRadians), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	for (int i = 2; i >= 0; i--)
	{
		Vector4f lTemp4VVector;
		int liIndex = liTriangleColNo - i;
		lTemp4VVector << V.col(liIndex)(0), V.col(liIndex)(1), 0, 1;
		Vector4f lFinalVectorValue = lRotationMatrix * lTemp4VVector;
		V.col(liIndex) << lFinalVectorValue(0), lFinalVectorValue(1);
	}

	//double lNewCentroidX = (V.col(liTriangleColNo - 2)(0) + V.col(liTriangleColNo - 1)(0) + V.col(liTriangleColNo)(0)) / 3.0;
	//double lNewCentroidY = (V.col(liTriangleColNo - 2)(1) + V.col(liTriangleColNo - 1)(1) + V.col(liTriangleColNo)(1)) / 3.0;

	//printCentroid(2);

	mdChangeTriangleX = lCentroidX - 0;
	mdChangeTriangleY = lCentroidY - 0;

	updateTranslationTriangle(miSelectedTriangle);

	//printCentroid(3);

	refreshMembers();

}


void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates
	double xcworld = ((xpos / double(width)) * 2) - 1;
	double ycworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw

	double ldScale = 1.0;
	if (mdZoom >= 0)
	{
		for (int z = 0; z < mdZoom; z++)
		{
			ldScale = ldScale - (ldScale * 0.2);
		}
	}
	else
	{
		for (int z = mdZoom; z < 0; z++)
		{
			ldScale = ldScale + (ldScale * 0.2);
		}
	}
	Matrix4f lInverseViewMatrix;
	lInverseViewMatrix << ldScale, 0, 0, -mdLR,
		0, ldScale, 0, -mdUD,
		0, 0, 1, 0,
		0, 0, 0, 1;

	Vector4f lTempVector;
	lTempVector << xcworld, ycworld, 0, 1;

	Vector4f lFinalVector = lInverseViewMatrix * lTempVector;

	double xworld = lFinalVector(0);
	double yworld = lFinalVector(1);

	//V.conservativeResize(V.rows(), mStepCounter + 1);
	if (miCurrentMode == 1)
	{
		V.col(miStepCounter) << xworld, yworld;
	}

	if (miCurrentMode == 2 && mbIsPressed == 1 && miSelectedTriangle != -1)
	{
		mdChangeTriangleX = xworld - mdSelectedTriangleX;
		mdChangeTriangleY = yworld - mdSelectedTriangleY;
		updateTranslationTriangle(miSelectedTriangle);
		mdSelectedTriangleX = xworld;
		mdSelectedTriangleY = yworld;
	}

	VBO.update(V);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the position of the mouse in the window
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates
	double xcworld = ((xpos / double(width)) * 2) - 1;
	double ycworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw


	double ldScale = 1.0;
	if (mdZoom >= 0)
	{
		for (int z = 0; z < mdZoom; z++)
		{
			ldScale = ldScale - (ldScale * 0.2);
		}
	}
	else
	{
		for (int z = mdZoom; z < 0; z++)
		{
			ldScale = ldScale + (ldScale * 0.2);
		}
	}
	Matrix4f lInverseViewMatrix;
	lInverseViewMatrix << ldScale, 0, 0, -mdLR,
		0, ldScale, 0, -mdUD,
		0, 0, 1, 0,
		0, 0, 0, 1;

	Vector4f lTempVector;
	lTempVector  << xcworld, ycworld, 0, 1;

	Vector4f lFinalVector = lInverseViewMatrix * lTempVector;

	double xworld = lFinalVector(0);
	double yworld = lFinalVector(1);

	//Kept Pressed
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		mbIsPressed = 1;

		if (miCurrentMode == 2)
		{
			int liOldTriangle = miSelectedTriangle;
			miSelectedTriangle = getTriangleNumber(xworld, yworld);
			selectTriangleColor(liOldTriangle, miSelectedTriangle);
			if (miSelectedTriangle != -1)
			{
				mdSelectedTriangleX = xworld;
				mdSelectedTriangleY = yworld;
			}
		}
	}

	//One Click (Press and Release)
	static int lsiOldState = GLFW_RELEASE;
	int liNewState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

	if (liNewState == GLFW_RELEASE && lsiOldState == GLFW_PRESS) {

		if (miCurrentMode == 1)
		{
			V.conservativeResize(V.rows(), V.cols() + 1);
			C.conservativeResize(C.rows(), C.cols() + 1);
			V.col(miStepCounter) << xworld, yworld;
			C.col(miStepCounter) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
			miStepCounter++;
			V.col(miStepCounter) << xworld, yworld;
			C.col(miStepCounter) << mLPinkColor(0), mLPinkColor(1), mLPinkColor(2);
			//	std::cout << V << "\n";
		}

		if (miCurrentMode == 3)
		{
			miSelectedTriangle = getTriangleNumber(xworld, yworld);
			//std::cout << liNumber << "\n";
			if (miSelectedTriangle != -1)
			{
				deleteTriangle(miSelectedTriangle);
			}
		}

		if (miCurrentMode == 4)
		{
			miClosestVertex = getClosestVertex(xworld, yworld);
			//std::cout << liNumber << "\n";
		}

		refreshMembers();
		mbIsPressed = 0;
	}
	lsiOldState = liNewState;

	// Upload the change to the GPU
	VBO.update(V);
	VBO_C.update(C);
}




void backupV()
{
	V_bkp.conservativeResize(V.rows(), V.cols());
	for (int i = 0; i < V.rows(); i++)
	{
		for (int j = 0; j < V.cols(); j++)
		{
			V_bkp(i, j) = V(i, j);
		}
	}
}

void restoreOriginalV()
{
	for (int i = 0; i < V_bkp.rows(); i++)
	{
		for (int j = 0; j < V_bkp.cols(); j++)
		{
			V(i, j) = V_bkp(i, j);
		}
	}
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
	string lStart = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?><svg width = \"200\" height = \"200\" xmlns = \"http://www.w3.org/2000/svg\" ><rect x = \"0\" y = \"0\" width = \"200\" height = \"200\" fill = \"rgb(76, 76, 76)\" />";
	string lTriangleTemplateStart = "<polygon points=\"";
	string lTTM = "\"";
	string lTriangleTemplateEnd = " stroke=\"rgb(255, 51, 153)\" stroke-width=\"2\" fill=\"rgb(255, 51, 153)\" />";

	string lend = "</svg>";
	string lfinalstring = lStart;

	for (int i = 0; i < V.cols() - 1; i = i + 3)
	{
		string lValues = "";

		int x = convertToSVGCoordinatesX(V.col(i)(0));
		int y = convertToSVGCoordinatesY(V.col(i)(1));
		lValues = lValues + to_string(x) + "," + to_string(y) + ",";

		x = convertToSVGCoordinatesX(V.col(i + 1)(0));
		y = convertToSVGCoordinatesY(V.col(i + 1)(1));
		lValues = lValues + to_string(x) + "," + to_string(y) + ",";

		x = convertToSVGCoordinatesX(V.col(i + 2)(0));
		y = convertToSVGCoordinatesY(V.col(i + 2)(1));
		lValues = lValues + to_string(x) + "," + to_string(y);

		string lfinaltriangle = lTriangleTemplateStart + lValues + lTTM + lTriangleTemplateEnd;
		lfinalstring = lfinalstring + lfinaltriangle;
	}

	lfinalstring = lfinalstring + lend;

	ofstream file("output.svg");
	file << lfinalstring;
	file.close();

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// Update the position of the first vertex if the keys 1,2, or 3 are pressed
	if (action == GLFW_PRESS)
	{
		if (miStepCounter % 3 != 0)
		{
			int liExtra = miStepCounter - ((miStepCounter / 3) * 3);
			for (int i = 0;i < liExtra;i++)
			{
				removeColumn(V, miStepCounter + 1);
				miStepCounter--;
			}
		}
		if (miCurrentMode == 5 && key != GLFW_KEY_R)
		{
			restoreOriginalV();
		}
		switch (key)
		{
		case GLFW_KEY_I:
			std::cout << "I is pressed" << "\n";
			miCurrentMode = 1;
			break;
		case GLFW_KEY_O:
			std::cout << "O is pressed" << "\n";
			miCurrentMode = 2;
			break;
		case GLFW_KEY_P:
			std::cout << "P is pressed" << "\n";
			miCurrentMode = 3;
			break;
		case GLFW_KEY_H:
			std::cout << "H is pressed" << "\n";
			miRotationMode = 1;
			updateRotationTriangle();
			miRotationMode = -1;
			break;
		case GLFW_KEY_J:
			std::cout << "J is pressed" << "\n";
			miRotationMode = 2;
			updateRotationTriangle();
			miRotationMode = -1;
			break;
		case GLFW_KEY_K:
			std::cout << "K is pressed" << "\n";
			miScalingMode = 1;
			updateScalingTriangle();
			miScalingMode = -1;
			break;
		case GLFW_KEY_L:
			std::cout << "L is pressed" << "\n";
			miScalingMode = 2;
			updateScalingTriangle();
			miScalingMode = -1;
			break;
		case GLFW_KEY_C:
			std::cout << "C is pressed" << "\n";
			miCurrentMode = 4;
			break;
		case GLFW_KEY_1:
			std::cout << "1 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mGreenColor);
			break;
		case GLFW_KEY_2:
			std::cout << "2 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mBlueColor);
			break;
		case GLFW_KEY_3:
			std::cout << "3 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mWhiteColor);
			break;
		case GLFW_KEY_4:
			std::cout << "4 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mOrangeColor);
			break;
		case GLFW_KEY_5:
			std::cout << "5 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mDGreenColor);
			break;
		case GLFW_KEY_6:
			std::cout << "6 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mYellowColor);
			break;
		case GLFW_KEY_7:
			std::cout << "7 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mBlackColor);
			break;
		case GLFW_KEY_8:
			std::cout << "8 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mRedColor);
			break;
		case GLFW_KEY_9:
			std::cout << "9 is pressed" << "\n";
			changeVertexColor(miClosestVertex, mLPinkColor);
			break;
		case GLFW_KEY_A:
			std::cout << "A is pressed" << "\n";
			mdLR = mdLR + 0.2;
			break;
		case GLFW_KEY_D:
			std::cout << "D is pressed" << "\n";
			mdLR = mdLR - 0.2;
			break;
		case GLFW_KEY_W:
			std::cout << "W is pressed" << "\n";
			mdUD = mdUD - 0.2;
			break;
		case GLFW_KEY_S:
			std::cout << "S is pressed" << "\n";
			mdUD = mdUD + 0.2;
			break;
		case GLFW_KEY_KP_ADD:
			std::cout << "+ is pressed" << "\n";
			mdZoom++;
			break;
		case GLFW_KEY_KP_SUBTRACT:
			std::cout << "- is pressed" << "\n";
			mdZoom--;
			break;
		case GLFW_KEY_R:
			std::cout << "R is pressed" << "\n";
			backupV();
			miCurrentMode = 5;
			break;
		case GLFW_KEY_T:
			std::cout << "T is pressed" << "\n";
			createSVG();
			break;
		default:
			break;
		}
		refreshMembers();
	}
	// Upload the change to the GPU
	VBO.update(V);
	VBO_C.update(C);
}


float mfOldTime[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
int miKFCounter[] = { 0,0,0,0,0,0,0,0,0,0 };

void updateKeyFrame(int pTriangleNo, double pKFCounter)
{
	Eigen::MatrixXf lDestination(2, 3);
	lDestination << 0, -0.5, 0.5,
		0.5, -0.5, -0.5;

	int liTriangleColNo = (pTriangleNo * 3) - 1;
	int liDestCounter = 0;
	for (int i = 2; i >= 0; i--)
	{
		Vector2f lOVector, lDVector, lFVector;
		int liIndex = liTriangleColNo - i;
		lOVector << V_bkp.col(liIndex)(0), V_bkp.col(liIndex)(1);
		lDVector << lDestination.col(liDestCounter)(0), lDestination.col(liDestCounter)(1);
		lFVector = ((1 - pKFCounter) * lOVector) + (pKFCounter * lDVector);
		V.col(liIndex) << lFVector(0), lFVector(1);
		liDestCounter++;
	}

		VBO.update(V);
	//	VBO_C.update(C);
}


void part1to5(void)
{
	GLFWwindow* window;

	// Initialize the library
	if (!glfwInit())
		return;

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
	window = glfwCreateWindow(640, 480, "Rupesh's World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		fprintf(stderr, "Create Window failed: \n");
		return;
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

	// V.resize(2,3);
	// V << 0,  0.5, -0.5, 0.5, -0.5, -0.5;
	VBO.update(V);
	VBO_C.update(C);

	// Initialize the OpenGL Program
	// A program controls the OpenGL pipeline and it must contains
	// at least a vertex shader and a fragment shader to be valid
	Program program;
	const GLchar* vertex_shader =
		"#version 150 core\n"
		"in vec2 position;"
		"in vec3 color;"
		"out vec3 f_color;"
		"uniform mat4 View;"
		"void main()"
		"{"
		"    gl_Position = View * vec4(position, 0.0, 1.0);"
		"    f_color = color;"
		"}";
	const GLchar* fragment_shader =
		"#version 150 core\n"
		"in vec3 f_color;"
		"out vec4 outColor;"
		"uniform vec3 triangleColor;"
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

	// Save the current time --- it will be used to dynamically change the triangle color
	auto t_start = std::chrono::high_resolution_clock::now();

	// Register the keyboard callback
	glfwSetKeyCallback(window, key_callback);

	// Register the mouse callback

	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);

	// Update viewport
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window))
	{
		// Bind your VAO (not necessary if you have only one)
		VAO.bind();

		// Bind your program
		program.bind();

		// Set the uniform value depending on the time difference
	   // auto t_now = std::chrono::high_resolution_clock::now();
	  // float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
	  // std::cout << time << "\n";
	  //  glUniform3f(program.uniform("triangleColor"), (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);

		// Clear the framebuffer
		glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		int liFullTriangleCounter = miStepCounter / 3;
		int liTriangleStartCounter = 0;
		int liTriangleNo = 1;
		//int liTempTriangleNo = 0;

		double ldScale = 1.0;
		if (mdZoom >= 0)
		{
			for (int z = 0; z < mdZoom; z++)
			{
				ldScale = ldScale + (ldScale * 0.2);
			}
		}
		else
		{
			for (int z = mdZoom; z < 0; z++)
			{
				ldScale = ldScale - (ldScale * 0.2);
			}
		}
		
		Matrix4f lViewMatrix;
		Matrix4f lDefaultViewMatrix;
		lDefaultViewMatrix << ldScale, 0, 0, mdLR,
			0, ldScale, 0, mdUD,
			0, 0, 1, 0,
			0, 0, 0, 1;

		//For Triangles
		while (liTriangleNo <= liFullTriangleCounter)
		{
			//std::cout << mSelectedTriangle << "\n";
			//if (liTriangleNo == miSelectedTriangle)
			//{
			//	glUniform3f(program.uniform("triangleColor"), 0.2f, 0.6f, 1.0f);
			//}
			//else
			//{
			//	glUniform3f(program.uniform("triangleColor"), 1.0f, 0.2f, 0.6f);
			//}

			//part5
			if (miCurrentMode == 5)
			{
				auto t_now = std::chrono::high_resolution_clock::now();
				float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
				if ((time - mfOldTime[liTriangleNo]) >= 1)
				{
					mfOldTime[liTriangleNo] = time;
					miKFCounter[liTriangleNo]++;
					if (miKFCounter[liTriangleNo] > 10)
					{
						miKFCounter[liTriangleNo] = 0;
					}
					//std::cout << liTriangleNo << "-------------" << miKFCounter[liTriangleNo] / 10.0f << "\n";
					updateKeyFrame(liTriangleNo, miKFCounter[liTriangleNo] / 10.0);
				}
			}

			glUniformMatrix4fv(program.uniform("View"), 1, GL_FALSE, lDefaultViewMatrix.data());
			glDrawArrays(GL_TRIANGLES, liTriangleStartCounter, 3);
			liTriangleStartCounter = liTriangleStartCounter + 3;
			liTriangleNo++;
		}

		//For Points and Lines
		if (miStepCounter % 3 != 0)
		{
			//glUniform3f(program.uniform("triangleColor"), 1.0f, 0.2f, 0.6f);
			glUniformMatrix4fv(program.uniform("View"), 1, GL_FALSE, lDefaultViewMatrix.data());

			glDrawArrays(GL_POINTS, liTriangleStartCounter, 1);
			glEnable(GL_PROGRAM_POINT_SIZE);
			glPointSize(2.0);

			glLineWidth(2.0);
			glDrawArrays(GL_LINES, liTriangleStartCounter, 2);
		}

		if (miStepCounter % 3 == 2)
		{
			//glUniform3f(program.uniform("triangleColor"), 1.0f, 0.2f, 0.6f);
			glUniformMatrix4fv(program.uniform("View"), 1, GL_FALSE, lDefaultViewMatrix.data());

			glDrawArrays(GL_TRIANGLES, liTriangleStartCounter, 3);
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
	VBO_C.free();

	// Deallocate glfw internals
	glfwTerminate();
	//return 0;
}
//////-----------------------------------------------------part1to5--------------------------------//////

// VertexBufferObject wrapper
VertexBufferObject VBOBC;

// Contains the vertex positions
Eigen::MatrixXf VBC(2, 15);

int miStepCounterBC = 0;
int miCurrentModeBC = -1;
int miCurrentBPoint = -1;
boolean lbdoChange = 0;


double bezier(double pA,  double pB,  double pC,  double pD,  double pT)  
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

void updateCurve()
{
	for (int i = 0;i <= 10;i++)
	{
		double x = bezier(VBC.col(0)(0), VBC.col(1)(0), VBC.col(2)(0), VBC.col(3)(0), i / 10.0);
		double y = bezier(VBC.col(0)(1), VBC.col(1)(1), VBC.col(2)(1), VBC.col(3)(1), i / 10.0);
	//	std::cout << x << " " << y << "\n";
		VBC.col(i + 4) << x, y;
	}
}

void cursor_position_callback6(GLFWwindow* window, double xpos, double ypos)
{
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates
	double xworld = ((xpos / double(width)) * 2) - 1;
	double yworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw

	//V.conservativeResize(V.rows(), mStepCounter + 1);
	if (lbdoChange == 1 && miCurrentBPoint >= 0 && miCurrentBPoint <= 3)
	{
		VBC.col(miCurrentBPoint) << xworld, yworld;
		updateCurve();
	}

	VBOBC.update(VBC);
}

void mouse_button_callback6(GLFWwindow* window, int button, int action, int mods)
{
	// Get the position of the mouse in the window
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	// Get the size of the window
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	// Convert screen position to world coordinates
	double xworld = ((xpos / double(width)) * 2) - 1;
	double yworld = (((height - 1 - ypos) / double(height)) * 2) - 1; // NOTE: y axis is flipped in glfw

	// Update the position of the first vertex if the left button is pressed
	//if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	//    V.col(0) << xworld, yworld;

	static int lsiOldState = GLFW_RELEASE;
	int liNewState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

	//std::cout << liNewState;
	if (liNewState == GLFW_RELEASE && lsiOldState == GLFW_PRESS)
	{
		//std::cout << VBC;
		if (miStepCounterBC >= 0 && miStepCounterBC <= 3)
		{
			for (int l = miStepCounterBC; l < 4; l++)
			{
				VBC.col(l) << xworld, yworld;
			}
			miStepCounterBC++;
		}
	}

	updateCurve();

	lsiOldState = liNewState;

	// Upload the change to the GPU
	VBOBC.update(VBC);
}

void key_callback6(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		// Update the position of the first vertex if the keys 1,2, or 3 are pressed
		switch (key)
		{
		case  GLFW_KEY_B:
			std::cout << "P is pressed" << "\n";
			miCurrentModeBC = 6;
			break;
		case  GLFW_KEY_X:
			std::cout << "X is pressed" << "\n";
			if (lbdoChange == 0)
				lbdoChange = 1;
			else
				lbdoChange = 0;
			break;
		case GLFW_KEY_N:
			std::cout << "N is pressed" << "\n";
			miCurrentBPoint++;
			if (miCurrentBPoint > 3)
				miCurrentBPoint = 0;
			break;
		case  GLFW_KEY_M:
			std::cout << "M is pressed" << "\n";
			miCurrentBPoint--;
			if (miCurrentBPoint < 0)
				miCurrentBPoint = 3;
			break;
		default:
			break;
		}
	}
	// Upload the change to the GPU
	VBOBC.update(VBC);
}


int part6(void)
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
	window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
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
	VBOBC.init();

	// V.resize(2,13);
	VBC << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	//V << 0, 0.5, -1, 0, 1, 0;
	VBOBC.update(VBC);

	// Initialize the OpenGL Program
	// A program controls the OpenGL pipeline and it must contains
	// at least a vertex shader and a fragment shader to be valid
	Program program;
	const GLchar* vertex_shader =
		"#version 150 core\n"
		"in vec2 position;"
		"void main()"
		"{"
		"    gl_Position = vec4(position, 0.0, 1.0);"
		"}";
	const GLchar* fragment_shader =
		"#version 150 core\n"
		"out vec4 outColor;"
		"uniform vec3 triangleColor;"
		"void main()"
		"{"
		"    outColor = vec4(triangleColor, 1.0);"
		"}";

	// Compile the two shaders and upload the binary to the GPU
	// Note that we have to explicitly specify that the output "slot" called outColor
	// is the one that we want in the fragment buffer (and thus on screen)
	program.init(vertex_shader, fragment_shader, "outColor");
	program.bind();

	// The vertex shader wants the position of the vertices as an input.
	// The following line connects the VBO we defined above with the position "slot"
	// in the vertex shader
	program.bindVertexAttribArray("position", VBOBC);

	// Save the current time --- it will be used to dynamically change the triangle color
	auto t_start = std::chrono::high_resolution_clock::now();

	// Register the keyboard callback
	glfwSetKeyCallback(window, key_callback6);

	// Register the mouse callback
	glfwSetCursorPosCallback(window, cursor_position_callback6);
	glfwSetMouseButtonCallback(window, mouse_button_callback6);

	// Update viewport
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	// Loop until the user closes the window
	while (!glfwWindowShouldClose(window))
	{
		// Bind your VAO (not necessary if you have only one)
		VAO.bind();

		// Bind your program
		program.bind();

		// Set the uniform value depending on the time difference
		//auto t_now = std::chrono::high_resolution_clock::now();
		//float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
		//glUniform3f(program.uniform("triangleColor"), (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);
		glUniform3f(program.uniform("triangleColor"), 1.0f, 0.2f, 0.6f);

		// Clear the framebuffer
		glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);

		// Draw a triangle
	  //  glDrawArrays(GL_TRIANGLES, 0, 3);


		glDrawArrays(GL_POINTS, 0, 4);
		glEnable(GL_PROGRAM_POINT_SIZE);
		glPointSize(7.0);

		glLineWidth(7.0);
		glDrawArrays(GL_LINE_STRIP, 4, 11);

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();
	}

	// Deallocate opengl memory
	program.free();
	VAO.free();
	VBOBC.free();

	// Deallocate glfw internals
	glfwTerminate();
	return 0;
}


int main(void)
{
	part1to5();
	part6();
	return 0;
}