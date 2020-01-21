#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"functions.h"

int main()
{	
	double alpha_d, gamma_d, current_d, ms_d, eta_d, beta_d, vol_d;
	double sim_time, temp, mz_tol, mz_der_tol, time_step_d;
	double  hext_d[3], Mf_d[3], han_d[3], hde_d[3];
	double theta_d, phi_d;
	int no_ensembles, choice;
	
	FILE *fp = fopen("Constants.d","r");
	if(fp == NULL){
		printf("Get the file containing the Constants, named as -Constants.d- !!");	
		return EXIT_FAILURE;	
	} 
	fscanf(fp,"%*s %lf",&alpha_d);
		fscanf(fp,"%*s %lf",&gamma_d);
		fscanf(fp,"%*s %lf",&current_d);
		fscanf(fp,"%*s %lf",&ms_d);
		fscanf(fp,"%*s %lf",&eta_d);
		fscanf(fp,"%*s %lf",&beta_d);
		fscanf(fp,"%*s %lf",&vol_d);
		fscanf(fp,"%*s %lf",&hext_d[0]);
		fscanf(fp,"%*s %lf",&hext_d[1]);
		fscanf(fp,"%*s %lf",&hext_d[2]);
		fscanf(fp,"%*s %lf",&han_d[0]);
		fscanf(fp,"%*s %lf",&han_d[1]);
		fscanf(fp,"%*s %lf",&han_d[2]);
		fscanf(fp,"%*s %lf",&hde_d[0]);
		fscanf(fp,"%*s %lf",&hde_d[1]);
		fscanf(fp,"%*s %lf",&hde_d[2]);
		fscanf(fp,"%*s %lf",&Mf_d[0]);
		fscanf(fp,"%*s %lf",&Mf_d[1]);
		fscanf(fp,"%*s %lf",&Mf_d[2]);
		fscanf(fp,"%*s %lf",&theta_d);
		fscanf(fp,"%*s %lf",&phi_d);		
		fscanf(fp,"%*s %d",&no_ensembles);
		fscanf(fp,"%*s %lf",&sim_time);
		fscanf(fp,"%*s %d",&choice);
		fscanf(fp,"%*s %lf",&temp);
		fscanf(fp,"%*s %lf",&mz_tol);
		fscanf(fp,"%*s %lf",&mz_der_tol);
		fscanf(fp,"%*s %lf",&time_step_d);
	fclose(fp);
	
	/*
	FILE *fp = fopen("values.d", "r");
	fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf",&alpha_d, &gamma_d, &current_d, &ms_d, &eta_d, &beta_d, &vol_d);
	fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", hext_d, hext_d+1, hext_d+2,han_d, han_d+1, han_d+2,hde_d, hde_d+1, hde_d+2, Mf_d, Mf_d+1, Mf_d+2);
	fscanf(fp, "%lf%lf", &theta_d, &phi_d);
	*/
	dostuff(alpha_d, gamma_d, current_d, ms_d, eta_d, beta_d, vol_d, hext_d, han_d, hde_d, Mf_d, theta_d, phi_d, no_ensembles,sim_time,choice, temp, mz_tol, mz_der_tol, time_step_d );
}
