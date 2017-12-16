#include "FUNC_CLAS.h"

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
	case 'S':
		base1 = VEC3D(1, 0, 0); base2 = VEC3D(0, 1, 0); base3 = VEC3D(0, 0, 1);
		break;
	case 'B':
		base1 = VEC3D(2, 0, 0); base2 = VEC3D(0, 2, 0); base3 = VEC3D(1, 1, 1);
		break;
	case 'F':
		base1 = VEC3D(1, 1, 0); base2 = VEC3D(1, 0, 1); base3 = VEC3D(0, 1, 1);
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


double VpsuFourier(VEC3D *q1, VEC3D *q2)
{
	VEC3D qq = *q1 - *q2;
	return -2.0*cos(qq.Norm() * 2)*exp(-qq.NormSqr()/5)/(qq.NormSqr() +.001);
}

#define AXISSEP 20
void DrawAxis(HDC hdc)
{
	MoveToEx(hdc, WIDTH + HBASE, BASE, 0);
	LineTo(hdc, HBASE, BASE);
	MoveToEx(hdc, HBASE, BASE + BASEADD, 0);
	LineTo(hdc, HBASE, 0);
	for (int i = HBASE; i < WIDTH + HBASE; i += 3)
		for (int j = 0; j < BASE + BASEADD; j += AXISSEP)
			SetPixel(hdc, i, j, RGB(128, 0, 255));
}

bool HDCToFile(const char* FilePath, HDC Context, RECT Area)
{
	uint16_t BitsPerPixel = 24;
	uint32_t Width = Area.right - Area.left;
	uint32_t Height = Area.bottom - Area.top;
	BITMAPINFO Info;
	BITMAPFILEHEADER Header;
	memset(&Info, 0, sizeof(Info));
	memset(&Header, 0, sizeof(Header));
	Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	Info.bmiHeader.biWidth = Width;
	Info.bmiHeader.biHeight = Height;
	Info.bmiHeader.biPlanes = 1;
	Info.bmiHeader.biBitCount = BitsPerPixel;
	Info.bmiHeader.biCompression = BI_RGB;
	Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3);
	Header.bfType = 0x4D42;
	Header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	char* Pixels = NULL;
	HDC MemDC = CreateCompatibleDC(Context);
	HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
	DeleteObject(SelectObject(MemDC, Section));
	BitBlt(MemDC, 0, 0, Width, Height, Context, Area.left, Area.top, SRCCOPY);
	DeleteDC(MemDC);
	FILE *hFile = fopen(FilePath, "wb");
	if (hFile)
	{
		fwrite((char*)&Header, sizeof(char), sizeof(Header) / sizeof(char), hFile);
		fwrite((char*)&Info.bmiHeader, sizeof(char), sizeof(Info.bmiHeader) / sizeof(char), hFile);
		fwrite(Pixels, sizeof(char), (((BitsPerPixel * Width + 31) & ~31) / 8) * Height / sizeof(char), hFile);
		fclose(hFile);
		DeleteObject(Section);
		return true;
	}
	DeleteObject(Section);
	return false;
}
