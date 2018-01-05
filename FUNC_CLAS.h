#ifndef FUNC_MAIN_H_H_H_H
#define FUNC_MAIN_H_H_H_H

#define _CRT_SECURE_NO_WARNINGS
#define HBASE 55
#define WIDTH 720
#define BASE 520
#define SCALE 500
#define BASEADD 230

#include <iostream>
//#include <iomanip>
//#include <complex>
#include <vector>
#include <cmath>
#include <Windows.h>
#include <mkl.h>

//using namespace std::complex_literals;
using namespace std;
//const std::complex<double> I = 1i;
const double pi = 3.1415926535897932385;
const double bohradiuspm = 52.917721092;
extern double lc;
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
void DrawAxis(HDC hdc);
void InitVariables();
double DivXBesselIN(int l, double alph, double qk);
double MulXBesselIN(int l, double alph, double qk);
double DivXBesselIL(int l, double alph, double qk);
double MulXBesselIL(int l, double alph, double qk);
double DivXBesselI(int l, double alph, double qk);
double MulXBesselI(int l, double alph, double qk);
double QI(double alph, double qq, double kk, int l);
double RI(double alph, double qq, double kk, int l);
double MXE(int ll, double qq, double kk);
int VpsuFourierNC(vector<VEC3D> *rcp, VEC3D *k, double *hmx, int n);
#endif

