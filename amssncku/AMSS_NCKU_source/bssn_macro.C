#include "bssn_macro.h"
#include <iostream>
#include <fstream>
#include <cstring>
using namespace std;

int compare_two_file(char *fname1, char *fname2, int data_num)
{
	// read file
	fstream file1(fname1, ios_base::in);
	fstream file2(fname2, ios_base::in);
	double *d1, *d2;
	d1 = (double *)malloc(sizeof(double) * data_num);
	d2 = (double *)malloc(sizeof(double) * data_num);

	for (int i = 0; i < data_num; ++i)
	{
		file1.read((char *)(d1 + i), sizeof(double));
		file2.read((char *)(d2 + i), sizeof(double));
	}

	// compare data
	bool is_match = true;
	for (int i = 0; i < data_num; ++i)
	{
		if (d1[i] != d2[i])
		{
			is_match = false;
			cout << "miss match at position " << i << endl;
			break;
		}
	}
	if (is_match)
		cout << "Result is right." << endl;

	free(d1);
	free(d2);
	file1.close();
	file2.close();
	return 0;
}
void printMatrix(int ftag1, int ftag2, double *d1, double *d2, int ord)
{
	char fname1[32];
	char fname2[32];
	// char ftag1[32]; char ftag2[32];
	// sprintf(ftag1,"%d",ftag1);
	strcpy(fname1, "matrix_f.show");
	// strcat(fname1,ftag1);

	// sprintf(ftag2,"%d",ftag2);
	strcpy(fname2, "matrix_g.show");
	// strcat(fname2,ftag2);

	ofstream fout0, fout1, fout2;
	fout1.open(fname1);
	fout2.open(fname2);

	for (int k = 0; k < 65; k++)
	{
		fout1 << "---------square " << k << " ----------" << endl;
		fout2 << "---------square " << k << " ----------" << endl;
		for (int j = 0; j < 67 + ord * 2; j++)
		{
			for (int i = 0; i < 67 + ord * 2; i++)
			{
				fout1 << d1[i + j * (67 + ord * 2) + k * ((67 + ord * 2) * (67 + ord * 2))] << ' ';
				fout2 << d2[i + j * (67 + ord * 2) + k * ((67 + ord * 2) * (67 + ord * 2))] << ' ';
				// fout1<<test_output_g[i+j*(cg->shape[0]) + k*(_2d_size)] <<' ';
				// fout2<<test_fh_f    [i+j*(cg->shape[0]) + k*(_2d_size)] <<' ';
			}
			fout1 << endl;
			fout2 << endl;
		}
	}
}

int compare_result(int ftag1, double *d2, int data_num)
{
	// read file
	char fname1[32];
	char ftag[32];
	// itoa(filetag,ftag,10);
	sprintf(ftag, "%d", ftag1);
	strcpy(fname1, "matrix_f.out");
	strcat(fname1, ftag);

	fstream file1(fname1, ios_base::in);
	double *d1;
	d1 = (double *)malloc(sizeof(double) * data_num);

	for (int i = 0; i < data_num; ++i)
	{
		file1.read((char *)(d1 + i), sizeof(double));
	}

	// compare data
	bool is_match = true;
	double delta;
	for (int i = 0; i < data_num; ++i)
	{
		delta = d1[i] - d2[i];
		if (delta < 0)
			delta = -delta;
		if (delta > 1e-14)
		{
			is_match = false;
			cout << fname1 << "::miss match at position " << i << endl;
			break;
		}
		// if(i<100 && i>50)
		//	cout<<d1[i]<<" "<<d2[i]<<endl;
	}
	if (is_match)
		cout << ftag1 << "::matched." << endl;

	if (ftag1 == 0)
	{
		printMatrix(1, 2, d1, d2, 3);
	}
	free(d1);
	file1.close();
	return 0;
}
