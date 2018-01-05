#include "FUNC_CLAS.h"

const int LMX = 3;

double Omega_c,//Volumn of the cell
lc,
Z_V,
c_core_1, c_core_2,
alpha_core_1, alpha_core_2,
alpha[LMX][6] =
{ { 1.92, 2.10, 2.39 }, { 0.82, 1.13, 1.51 }, { 1.36, 1.59, 1.77 } }
,
A[LMX][6] =
{ { -2.6670, .7075, .0251, .0608, -.0134, -.0161 },
{ -1.5706, -.2352, .0327, .0262, -.0090, -.0047 },
{ -.2574, -.5358, -.0668, -.1835, .0187, .0551 } }
, vecscale;

VEC3D::VEC3D()
{
}

VEC3D::VEC3D(double xx, double yy, double zz)
{
	sx = xx; sy = yy; sz = zz;
}

VEC3D::~VEC3D()
{
}

VEC3D operator+(VEC3D v1, VEC3D v2)
{
	return VEC3D(v1.sx + v2.sx, v1.sy + v2.sy, v1.sz + v2.sz);
}

VEC3D operator-(VEC3D v1, VEC3D v2)
{
	return VEC3D(v1.sx - v2.sx, v1.sy - v2.sy, v1.sz - v2.sz);
}

VEC3D operator*(VEC3D v1, VEC3D v2)
{
	return VEC3D((v1.sy)*(v2.sz) - (v1.sz)*(v2.sy),
		(v1.sz)*(v2.sx) - (v1.sx)*(v2.sz),
		(v1.sx)*(v2.sy) - (v1.sy)*(v2.sx));
}

VEC3D operator*(double lbd, VEC3D v1)
{
	return VEC3D((v1.sx)*lbd,
		(v1.sy)*lbd,
		(v1.sz)*lbd);
}

VEC3D operator*(VEC3D v1, double lbd)
{
	return VEC3D((v1.sx)*lbd,
		(v1.sy)*lbd,
		(v1.sz)*lbd);
}

double operator|(VEC3D v1, VEC3D v2)
{
	return (v1.sx)*(v2.sx)
		+ (v1.sy)*(v2.sy)
		+ (v1.sz)*(v2.sz);
}

double VEC3D::Norm()
{
	return sqrt(sx*sx + sy*sy + sz*sz);
}

double VEC3D::NormSqr()
{
	return sx*sx + sy*sy + sz*sz;
}

void VEC3D::Print()
{
	cout << '(' << sx << ',' << sy << ',' << sz << ')' << endl;
}

void GetWaveVector(char rltype, double cutOffErg, vector<VEC3D> *arry, int *n)
{
	int kmax = 2 * (1 + (int)sqrt(cutOffErg)), cont = 0, nsz = kmax*kmax*kmax;
	arry->resize(nsz);
	VEC3D base1, base2, base3, vec;
	switch (rltype)
	{
	case 'B':
		base1 = VEC3D(1, 1, 0)*(2 * pi / lc); base2 = VEC3D(0, 1, 1)*(2 * pi / lc); base3 = VEC3D(1, 0, 1)*(2 * pi / lc);
		break;
	case 'F':
		base1 = VEC3D(1, 1, -1)*(2 * pi / lc); base2 = VEC3D(1, -1, 1)*(2 * pi / lc); base3 = VEC3D(-1, 1, 1)*(2 * pi / lc);
		break;
	default:
		*n = 0;
		return;
	}
	for (int i = -kmax; i < kmax; i++)
		for (int j = -kmax; j < kmax; j++)
			for (int k = -kmax; k < kmax; k++)
			{
				vec = ((double)i*base1) + ((double)j*base2) + ((double)k*base3);
				if (vec.NormSqr() < cutOffErg)
				{
					(*arry)[cont++] = vec;
					if (cont == nsz)
					{
						nsz += 20;
						arry->resize(nsz);
					}
				}
			}
	arry->resize(nsz);
	*n = cont;
}

void InitVariables()
{
	double Qij[36], Sij[36], sqrtpi = sqrt(pi);

	lc = 404.95 / bohradiuspm;
	Omega_c = lc*lc*lc / 4.0;
	Z_V = 3.0;
	alpha_core_1 = 1.77; alpha_core_2 = .70;
	c_core_1 = 1.7905; c_core_2 = 1.0 - c_core_1;

	for (int l = 0; l < LMX; l++)
	{
		//Duplicate C
		for (int i = 0; i < 3; i++)
			alpha[l][i + 3] = alpha[l][i];

		//Generate Sij
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (i < 3 && j < 3)	Sij[i * 6 + j] = sqrtpi / (4.0*pow(alpha[l][i] + alpha[l][j], 1.5));
				else if (i > 2 && j > 2) Sij[i * 6 + j] = 0.9375*sqrtpi / pow(alpha[l][i] + alpha[l][j], 3.5);
				else Sij[i * 6 + j] = 0.375*sqrtpi / pow(alpha[l][i] + alpha[l][j], 2.5);
			}
		}

		//Generate Qij
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				if (i > j) Qij[i * 6 + j] = 0;
				else if (i == j)
				{
					Qij[i * 6 + j] = Sij[i * 6 + j];
					for (int kk = 0; kk < i; kk++)
						Qij[i * 6 + j] -= (Qij[kk * 6 + i] * Qij[kk * 6 + i]);
					if (Qij[i * 6 + j] < 0) cout << "Error in Qij!\n";
					else Qij[i * 6 + j] = sqrt(Qij[i * 6 + j]);
				}
				else
				{
					Qij[i * 6 + j] = Sij[i * 6 + j];
					for (int kk = 0; kk < i; kk++)
						Qij[i * 6 + j] -= (Qij[kk * 6 + i] * Qij[kk * 6 + j]);
					Qij[i * 6 + j] /= Qij[i * 6 + i];
				}
			}
		}

		//Get Inverse;
		cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 6, Qij, 6, A[l], 1);
		for (int i = 0; i < 6; i++) A[l][i] = -A[l][i];
	}
}

double LegendreP(int n, double x)
{
	switch (n)
	{
	case 0:
		return 1;
		break;
	case 1:
		return x;
		break;
	case 2:
		return 1.5*x*x - .5;
		break;
	case 3:
		return x*(2.5*x*x - 1.5);
		break;
	default:
		cout << "LegendreP Error!";
		return 0;
		break;
	}
	return 0;
}

double DivXBesselIN(int l, double alph, double qk)
//pi/sqrt(qk)I(l+1/2,qk/(2 alph))
{
	double qko2a = qk / 2 / alph;
	switch (l)
	{
	case 0:
		return 2 * sqrt(alph*pi)*sinh(qko2a) / qk;
		break;
	case 1:
		return 2 * sqrt(alph*pi)*(qk*cosh(qko2a) - 2 * alph*sinh(qko2a)) / qk / qk;
		break;
	case 2:
		return 2 * sqrt(alph*pi)*(-6 * qk*alph*cosh(qko2a) + (qk*qk + 12 * alph*alph)*sinh(qko2a)) / qk / qk / qk;
		break;
	case 3:
		return 2 * sqrt(alph*pi)*(qk *(qk*qk + 60 * alph*alph)*cosh(qko2a) - 12 * alph*(qk*qk + 10 * alph*alph)*sinh(qko2a)) / qk / qk / qk / qk;
		break;
	default:
		break;
	}
	return 0;
}

double MulXBesselIN(int l, double alph, double qk)
//pi*sqrt(qk)*I(l-1/2,qk/(2 alph))
{
	double qko2a = qk / 2 / alph;
	switch (l)
	{
	case 0:
		return  2 * sqrt(alph*pi)*cosh(qko2a);
		break;
	case 1:
		return 2 * sqrt(alph*pi)*sinh(qko2a);
		break;
	case 2:
		return 2 * sqrt(alph*pi)*(qk*cosh(qko2a) - 2 * alph*sinh(qko2a)) / qk;
		break;
	case 3:
		return 2 * sqrt(alph*pi)*(-6 * qk*alph*cosh(qko2a) + (qk*qk + 12 * alph*alph)*sinh(qko2a)) / qk / qk;
		break;
	default:
		break;
	}
	return 0;
}

double DivXBesselIL(int l, double alph, double qk)
//pi/sqrt(qk)I(l+1/2,qk/(2 alph))
{
	double temp;
	switch (l)
	{
	case 0:
		temp = pow(qk / alph, 2);
		return sqrt(pi / alph)*(1 + temp / 24 + temp*temp / 1920);
		break;
	case 1:
		temp = qk / alph;
		return sqrt(pi / alph)*(temp / 6 + pow(temp, 3) / 240);
		break;
	case 2:
		temp = pow(qk / alph, 2);
		return sqrt(pi / alph)*(temp / 60 + temp*temp / 3360);
		break;
	case 3:
		temp = qk / alph;
		return sqrt(pi / alph)*pow(temp, 3) / 840;
		break;
	default:
		break;
	}
	return 0;
}

double MulXBesselIL(int l, double alph, double qk)
//pi*sqrt(qk)*I(l-1/2,qk/(2 alph))
{
	double temp;
	switch (l)
	{
	case 0:
		temp = pow(qk / alph, 2);
		return sqrt(pi * alph)*(2 + temp / 4 + temp*temp / 192);
		break;
	case 1:
		temp = qk / alph;
		return sqrt(pi * alph)*(temp + temp*temp*temp / 24);
		break;
	case 2:
		temp = pow(qk / alph, 2);
		return sqrt(pi * alph)*(temp / 6 + temp*temp / 240);
		break;
	case 3:
		temp = qk / alph;
		return sqrt(pi * alph)*temp*temp*temp / 60;
		break;
	default:
		break;
	}
	return 0;
}

double DivXBesselI(int l, double alph, double qk)
{
	if (qk / alph > .08) return DivXBesselIN(l, alph, qk);
	else return DivXBesselIL(l, alph, qk);
	return 0;
}

double MulXBesselI(int l, double alph, double qk)
{
	if (qk / alph > .08) return MulXBesselIN(l, alph, qk);
	else return MulXBesselIL(l, alph, qk);
	return 0;
}

double QI(double alph, double qq, double kk, int l)
{
	return exp(-(qq*qq + kk*kk) / 4 / alph)
		*DivXBesselI(l, alph, qq*kk) / 4 / alph;
}

double RI(double alph, double qq, double kk, int l)
{
	return exp(-(qq*qq + kk*kk) / 4 / alph) / 4 / alph / alph
		*(MulXBesselI(l, alph, qq*kk) / 2 / alph
			- ((qq*qq + kk*kk) / 4 / alph + l - .5)*DivXBesselI(l, alph, qq*kk));
}

double MXE(int ll, double qq, double kk)
{
	double xxx = 0;
	for (int ig = 0; ig < 3; ig++)
	{
		xxx += A[ll][ig] * QI(alpha[ll][ig], qq, kk, ll);
		xxx += A[ll][ig + 3] * RI(alpha[ll][ig], qq, kk, ll);
	}
	return xxx;
}

int VpsuFourierNC(vector<VEC3D> *rcp, VEC3D *k, double *hmx, int n)
{
	VEC3D q;
	//VEC3D dtau(lc / 8, lc / 8, lc / 8);
	double qq, temp;
	//Local
	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
		{
			qq = ((*rcp)[j] - (*rcp)[i]).NormSqr();
			temp = -4.0* pi*Z_V / (qq*Omega_c)*(c_core_1*exp(-qq / (4.0*alpha_core_1)) + c_core_2*exp(-qq / (4.0*alpha_core_2)));
			hmx[i + j*n] = temp;// *2 * cos(dtau | ((*rcp)[j] - (*rcp)[i]));
		}

	for (int i = 0; i < n; i++)
		hmx[i + i*n] = 0;

	//Non local
	for (int li = 0; li < LMX; li++)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++)
			{
				temp = 4.0 * pi / Omega_c*(2 * li + 1)*MXE(li, ((*rcp)[i] + *k).Norm(), ((*rcp)[j] + *k).Norm());
				temp *= LegendreP(li, (((*rcp)[i] + *k) | ((*rcp)[j] + *k)) / ((((*rcp)[j] + *k).Norm())*(((*rcp)[i] + *k).Norm()) + 1e-20));
				hmx[i + j*n] += temp;// *2 * cos(dtau | ((*rcp)[j] - (*rcp)[i]));
				//	if (temp != temp) system("pause");
			}
	}/**/


	//Kinetic
	for (int i = 0; i < n; i++)
		hmx[i + i*n] += ((*rcp)[i] + *k).NormSqr() / 2;

	return 0;
}

#define AXISSEP 50
void DrawAxis(HDC hdc)
{
	char no[30];
	MoveToEx(hdc, WIDTH + HBASE, BASE, 0);
	LineTo(hdc, HBASE, BASE);
	MoveToEx(hdc, HBASE, BASE + BASEADD, 0);
	LineTo(hdc, HBASE, 0);
	for (int j = BASE; j > 0; j -= AXISSEP)
	{
		for (int i = HBASE; i < WIDTH + HBASE; i += 3)
			SetPixel(hdc, i, j, RGB(128, 0, 255));
	}
	for (int j = BASE; j > 0; j -= AXISSEP * 2)
		TextOutA(hdc, HBASE - 45, j - 8, no, sprintf(no, "%.2f", 2.0*(BASE - j) / SCALE));
}