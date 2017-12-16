#ifndef FUNC_MAIN_H_H_H_H
#define FUNC_MAIN_H_H_H_H

#define _CRT_SECURE_NO_WARNINGS
#define HBASE 30
#define WIDTH 720
#define BASE 520
#define SCALE 100
#define BASEADD 30


#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <cmath>
#include <Windows.h>
#include <mkl.h>

using namespace std::complex_literals;
using namespace std;
const std::complex<double> I = 1i;
const double pi = 3.1415926535897932385;

class VEC3D
{
public:
	VEC3D();
	VEC3D(double xx, double yy, double zz);
	~VEC3D();
	double Norm();
	double NormSqr();
	double sx, sy, sz;
	void Print();
};

VEC3D operator+(VEC3D v1, VEC3D v2);
VEC3D operator+(VEC3D v1, VEC3D v2);
VEC3D operator-(VEC3D v1, VEC3D v2);
VEC3D operator*(VEC3D v1, VEC3D v2);
VEC3D operator*(double lbd, VEC3D v1);
VEC3D operator*(VEC3D v1, double lbd);
double operator|(VEC3D v1, VEC3D v2);

void GetWaveVector(char rltype, double cutOffErg, vector<VEC3D> *arry, int *n);
double VpsuFourier(VEC3D *q1, VEC3D *q2);
bool HDCToFile(const char* FilePath, HDC Context, RECT Area);
void DrawAxis(HDC hdc);
#endif

