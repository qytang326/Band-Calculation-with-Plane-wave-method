#include "FUNC_CLAS.h"

HDC hdc;

void ErgBand(double cutErg, VEC3D *stt, int *pltpts, int nplt, wchar_t **name)
{
	const int nbd = 25;
	FILE *fp = fopen("eg.txt", "w");
	int n = 0;
	vector<VEC3D> rcplts;
	double *hmtn;
	double *diag, *edg, *tau;
	VEC3D dk, ks(0.001, 0.001, 0.001);
	double dx, xx = HBASE;
	for (int i = 0; i < nplt; i++)
		n += pltpts[i];
	dx = (double)WIDTH / n;

	GetWaveVector('F', cutErg, &rcplts, &n);

	cout << "N=" << n << endl;
	hmtn = new double[n*n];
	diag = new double[n];
	edg = new double[n];
	tau = new double[n];

	RECT rc = { 0, 0, 2 * HBASE + WIDTH, BASE + BASEADD };
	FillRect(hdc, &rc, (HBRUSH)(COLOR_WINDOW + 1));
	VpsuFourierNC(&rcplts, &ks, hmtn, n);

	LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', n, hmtn, n, diag, edg, tau);
	while (LAPACKE_dsterf(n, diag, edg));

	double offset = diag[0];

	for (int j = 0; j < nbd; j++)
	{
		SetPixel(hdc, (int)xx, BASE - (int)(SCALE*(diag[j] - offset)), 0);
		fprintf(fp, "%f,", diag[j]);
	}
	fprintf(fp, "%f\n", diag[nbd]);

	DrawAxis(hdc);

	for (int np = 0; np < nplt; np++)
	{
		dk = (stt[np + 1] - stt[np])*(1.0 / (pltpts[np]));
		TextOutW(hdc, (int)xx + 10, BASE + 10, name[np], lstrlenW(name[np]));
		for (int i = 0; i < pltpts[np]; i++)
		{
			fprintf(fp, "0");
			ks = (i + 1)*dk + stt[np];

			VpsuFourierNC(&rcplts, &ks, hmtn, n);

			LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', n, hmtn, n, diag, edg, tau);
			while (LAPACKE_dsterf(n, diag, edg));

			//if ((diag[0] - offset) < 0) system("pause");
			for (int j = 0; j < nbd; j++)
			{
				SetPixel(hdc, (int)xx, BASE - (int)(SCALE*(diag[j] - offset)), 0);
				fprintf(fp, "%f,", diag[j]);
			}
			fprintf(fp, "%f\n", diag[nbd]);
			xx += dx;

		}
		MoveToEx(hdc, (int)(xx - dx), 0, 0);
		LineTo(hdc, (int)(xx - dx), BASE + BASEADD);
	}
	TextOutW(hdc, (int)xx + 10, BASE + 10, name[nplt], lstrlenW(name[nplt]));

	delete[] hmtn;
	delete[] diag;
	delete[] edg;
	delete[] tau;
	fclose(fp);
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

	InitVariables();
	VEC3D vst[] = { VEC3D(0, 0, 0), VEC3D(0, 1, 0), VEC3D(.5, 1, 0), VEC3D(.5, .5, .5), VEC3D(0, 0, 0), VEC3D(.75, .75, 0) };
	for (int i = 0; i < 6; i++) vst[i] = vst[i] * (2 * pi / lc);
	wchar_t *name[] = { TEXTGAMMA, L"X", L"W", L"L", TEXTGAMMA, L"K" };
	//int pn[] = { 20, 10, 15, 20, 20 };
	int pn[] = { 200, 100, 150, 200, 200 };
	ErgBand(15, vst, pn, 5, name);
	return system("pause");
}
