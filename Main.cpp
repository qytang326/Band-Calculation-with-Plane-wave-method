#include "FUNC_CLAS.h"

HDC hdc;

void ErgBand(double cutErg, VEC3D *stt, int *pltpts, int nplt, wchar_t **name)
{
	int n = 0;
	vector<VEC3D> rcplts;
	double *hamiltonian, *hmtn;
	double *diag, *edg, *tau;
	VEC3D dk, ks;

	double dx, xx = HBASE;
	for (int i = 0; i < nplt; i++)
		n += pltpts[i];
	dx = (double)WIDTH / n;

	GetWaveVector('B', cutErg, &rcplts, &n);
	cout << "N=" << n << endl;
	hamiltonian = new double[n*n];
	hmtn = new double[n*n];
	diag = new double[n];
	edg = new double[n - 1];
	tau = new double[n];

	RECT rc = { 0,0,2 * HBASE + WIDTH,BASE + BASEADD };
	FillRect(hdc, &rc, (HBRUSH)(COLOR_WINDOW + 1));

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			hmtn[i + j*n] = hamiltonian[i + j*n] = VpsuFourier(&rcplts[i], &rcplts[j]);
	for (int i = 0; i < n; i++)
		hmtn[i + i*n] = rcplts[i].NormSqr();

	LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', n, hmtn, n, diag, edg, tau);
	while (LAPACKE_dsterf(n, diag, edg));
	double offset = diag[0];

	DrawAxis(hdc);

	for (int np = 0; np < nplt; np++)
	{
		dk = (stt[np + 1] - stt[np])*(1.0 / (pltpts[np] - 1));
		TextOutW(hdc, (int)xx + 10, BASE + 10, name[np], lstrlenW(name[np]));
		for (int i = 0; i < pltpts[np]; i++)
		{
			ks = i*dk + stt[np];
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					hmtn[i + j*n] = hamiltonian[i + j*n];
			for (int i = 0; i < n; i++)
				hmtn[i + i*n] = (ks + rcplts[i]).NormSqr();

			LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', n, hmtn, n, diag, edg, tau);
			while (LAPACKE_dsterf(n, diag, edg));

			for (int i = 0; i < n; i++)
				SetPixel(hdc, (int)xx, BASE - (int)(SCALE*(diag[i] - offset)), 0);
			xx += dx;
		}
		MoveToEx(hdc, (int)(xx - dx), 0, 0);
		LineTo(hdc, (int)(xx - dx), BASE + BASEADD);
	}
	TextOutW(hdc, (int)xx + 10, BASE + 10, name[nplt], lstrlenW(name[nplt]));
	delete[] hamiltonian;
	delete[] hmtn;
	delete[] diag;
	delete[] edg;
	delete[] tau;
	HDCToFile("ErgBand.bmp", hdc, rc);
}

INT_PTR CALLBACK WndProc(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	if (message == WM_INITDIALOG)
	{
		hdc = GetDC(hDlg);
		SelectObject(hdc, CreateFontW(0, 0, 0, 0, FW_BOLD, 0, 0, 0,
			GREEK_CHARSET, OUT_CHARACTER_PRECIS, CLIP_CHARACTER_PRECIS,
			DEFAULT_QUALITY, FF_DONTCARE, L"Times New Roman"));
		return (INT_PTR)1;
	}
	else if (message == WM_CLOSE)
	{
		ReleaseDC(hDlg, hdc);
		EndDialog(hDlg, 0);
	}
	return (INT_PTR)0;
}

void WINAPI OpenWindow()
{
	DialogBox(0, MAKEINTRESOURCE(101), 0, WndProc);
}

#define TEXTGAMMA L"\x393"

int main()
{
	CloseHandle(CreateThread(0, 0, (LPTHREAD_START_ROUTINE)OpenWindow, 0, 0, 0));
	Sleep(1000);
	//SC
	//VEC3D vst[] = { VEC3D(0,0,0), VEC3D(0,0,.5), VEC3D(0,.5,.5), VEC3D(.5,.5,.5), VEC3D(0,0,0) };
	//wchar_t *name[] = { TEXTGAMMA, L"M",L"X",L"R", TEXTGAMMA };
	//int pn[] = { 400, 200, 200, 400 };
	//ErgBand(14, vst, pn, 4, name);

	//BCC
	VEC3D vst[] = { VEC3D(0,0,0), VEC3D(0,1,0), VEC3D(.5,1,0), VEC3D(.5,.5,.5), VEC3D(0,0,0), VEC3D(.75,.75,0) };
	wchar_t *name[] = { TEXTGAMMA, L"X",L"W",L"L", TEXTGAMMA , L"K" };
	int pn[] = { 400, 200, 200, 400,400 };
	ErgBand(50, vst, pn, 5, name);
	return system("pause");
}
