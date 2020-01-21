#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define CLEARBUF() char ch; while((ch=getchar())!= EOF && ch != '\n');

void computeAvgVarSkew(char *, char *);

int main(int argc, char ** argv)
{
	char *readfile = malloc(300*sizeof(char));
	char *writefile = malloc(300*sizeof(char));
	readfile = argv[1];
	strcat(writefile, "avg_var_skew");
	strcat(writefile, readfile);	
	computeAvgVarSkew(readfile, writefile);
	printf(".....................Done!\n");

}

void computeAvgVarSkew(char *readfile, char *writefile)
{
	
	FILE *fp1 = fopen(readfile, "r");
	if(fp1==NULL)
	{
		printf("No such file\n");
		return;
	}
	FILE *fp2 = fopen(writefile, "w");
	double tsavg=0, tsqavg=0, tcubavg=0; int no;
	double ti;

	while(!feof(fp1))
	{
		fscanf(fp1, "%d%lf", &no, &ti);
		tsavg+=ti;
		tsqavg+=ti*ti;
		tcubavg+=ti*ti*ti;
	}
	tsavg = tsavg/no;
	tsqavg = tsqavg/no;
	tcubavg = tcubavg/no;
	double sigma_t = sqrt(tsqavg - tsavg*tsavg);

	double skewness_t = cbrt(tcubavg - 3*sigma_t*tsavg + 2*tsavg*tsavg*tsavg);

	fprintf(fp2, "%.4f\t%.4f\t%.4f", 3.483*tsavg, 3.483*sigma_t, 3.483*skewness_t);
	
	fclose(fp1);
	fclose(fp2);
	return;

}
