#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "macrodef.h"
extern "C"
{
#ifdef fortran1
	void writefile_f
#endif
#ifdef fortran2
		void WRITEFILE_F
#endif
#ifdef fortran3
		void
		writefile_f_
#endif
		(int &filetag, double *matrix, int &msize)
	{
		char fname[32];
		char ftag[32];
		// itoa(filetag,ftag,10);
		sprintf(ftag, "%d", filetag);
		strcpy(fname, "matrix_f.out");
		strcat(fname, ftag);

		/*printf("-------------called-------------");
		printf(fname);
		printf("\n");
		printf("int tag %d\n",filetag);
		printf("int msize %d\n",msize);
		printf(ftag);*/

		printf("int msize %d\n", msize);

		FILE *fp;
		fp = fopen(fname, "wb");
		// char buffer[1024];
		// buffer[1023]='\0';
		// int bsize;

		if (fp == NULL)
		{
			printf("Open file failed.");
			exit(0);
		}

		// msize =  sizeof(double) * msize;
		fwrite(matrix, sizeof(double), msize, fp);

		fclose(fp);
		// return 0;
	}
}
