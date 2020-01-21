#include<stdio.h>
#include<stdlib.h>
#include"functions.h"
#include"gauss.c"
double alpha, gama, current, ms, eta, beta, vol, K;
double mz_tol, mz_der_tol;
double heff[3], mf[3], hext[3], han[3], hde[3];
double theta, phi;
int NO_TIME_STEPS, NO_ENSEMBLES;
double KB_TIMES_TEMP;

double ti, time_step;

double aj_var, bj_var;


long seed = -10000;
/*
 * The following function is the initialization function.
 * It takes paramters from the main file and stores them in global variables declared in the "functions.h" file.
 * The "_d" suffix in variable names below indicates that the values of the variables come with dimensions. i.e, they are real, physical values
 * The function below converts the values into dimensionless quantities by dividing with appropriate constants.
*/


void dostuff( double alpha_d, double gamma_d,double current_d, double ms_d, double eta_d, double beta_d, double vol_d, double *hext_d, double *han_d, double *hde_d,  double *Mf_d, double theta_d, double phi_d, int no_ensembles, double sim_time_d, int choice, double temp_d, double mz_tol_d, double mz_der_tol_d, double time_step_d )
{
	//Magnitude of the effective magnetic field. Its usage is to normalize the fixed magnetic field vector, Mf
	double mag_mf=0;
	//Dummy variable used for loops	
	int i;
	
	//Operation index. Is used to take in input the operation the user wishes to do.	
	int op=choice;
	//Character which upon input prompt if filled with 'y' or 'n' depending upon whether or not to shutdown pc after simulation
	//Note the user needs to be root for the system to shutdown after completion of simulation
	char shutdown;
	//Function pointer to point to function decided by the "op" variable. Both "op" and "wtd" are assigned values further in the program
	void (*wtd)();
	
	//Alpha is just the damping constant.
	alpha = alpha_d;
	//Constantof proportionality in dm/dt = gama * (m x Heff)
	gama = gamma_d;
	//The current supplied to give rise to Spin Transfer Torque and Field Like Torque
	current = current_d;
	//The magnitude of the soft magnetic field
	ms = ms_d;
	//Used for calculation of aj- the term corresponding to STT
	eta = eta_d;
	//Ratio of aj and bj. ie bj = beta*aj
	beta = beta_d;
	//Volume of the Magnetic Tunneling Junction
	vol = vol_d;
	
	//Initia values of theta and phi
	theta = theta_d;
	phi = phi_d;
	
	//Initial value of time and the size of the time step.


	/*  RATIONALE OF THE TIME STEP SIZE
	    -------------------------------
		
	 * The time step taken in real time is in Pico Second (10^-12 s)
	 * If t is a time-step then (t*4*pi*ms*gama) is the dimensionless time step
	 * 4*pi*ms*gama is equal to 0.285 * 10 ^ 12 s.
	 * Therefore 1 picoS is equal to 0.285 in dimensionless time.
	 */

	ti=0; time_step=4*M_PI*ms*gama*time_step_d;
	printf("Step size is %lf\n", time_step);

	/*
		RATIONALE OF THE NO OF TIME STEPS
		---------------------------------
	
		* The sim_time is the value of the actual period of time that the program is simulating for the system.
		* The number of time steps therefore is the sim_time divided by the step size.
	*/
	NO_TIME_STEPS = (int)(sim_time_d/time_step_d);
	printf("No of time_steps  is %d\n",NO_TIME_STEPS );


	/*  Initialize the number of Ensembles*/
	NO_ENSEMBLES = no_ensembles;
	printf("No of ensembles  is %d\n",NO_ENSEMBLES );

	/*  Initialize the KB_TIMES_TEMP value*/
	KB_TIMES_TEMP = KB*temp_d;
	printf("Time_step: %e \n", time_step_d);

	/*Initialize the mz tolerance and the derivative z tolerance*/
	mz_tol = mz_tol_d;
	mz_der_tol = mz_der_tol_d;
	
	//Constants contributing to the Effective Field
        // 1. hext ------> The external applied magnetic field
	// 2. han  ------> Field arising from anisotropy. Dependent upon geometry. In this case,in z direction only.
	// 3. hde  ------> Field arising from the magnetization of the fixed magnet, Mf. Dependent on geomtery. In this case, in z direction only.
	for(i=0; i<3; i++)
	{
		hext[i] = hext_d[i];
		han[i] = han_d[i];
		hde[i] = hde_d[i];
	}

	//Making the quantities dimensionless by dividing with (4*pi*ms)
	for(i=0; i<3; i++)
	{
		hext[i] = hext[i]/(4*M_PI*ms);
		han[i] = han[i]/(4*M_PI*ms);
		hde[i] = hde[i]/(4*M_PI*ms);
	}

	//Diving Mf with its magnitude to make it dimensionless and a unit vector
	
	for(i=0; i<3; i++)
	mag_mf = mag_mf + Mf_d[i]*Mf_d[i];

	mag_mf = sqrt(mag_mf);

	for(i=0; i<3; i++)
	mf[i] = Mf_d[i]/mag_mf;
	
	//Display the options of the operations
	displayPrompt();	
		
	printf("value of op is %d\n", op);
	//The following switch block examines the value of "op" and assigns the address of the appropriate function to "wtd".
	switch(op)
	{
		case 1:
		wtd = &heffOnly;
		break;
		
		case 2:
		wtd = &heffDamp;
		break;
		
		case 3:
		wtd = &heffDampStt;
		break;

		case 4:
		wtd = &heffDampSttFlt;
		break;

		case 5:
		wtd = &heffDampFluc;
		break;

		case 6:
		wtd = &heffDampSttFluc;
		break;

		case 7:
		wtd = &heffDampSttFltFluc;
		break;
		
		case 8:
		wtd = &flucOnly;
		break;
		
	
		default:
		wtd = &heffDampSttFltFluc;

	}
	printf("\n\n\nValues of theta and phi: %lf, %lf\n\n\n ", theta, phi);
	/*printf("Do you wish to shut down after completion?[y/n]\n");
	scanf("%c", &shutdown);*/
	

	clock_t tic = clock();
	(*wtd)();
	clock_t toc = clock();
	printf("Done computing\n");
	 
   	FILE *tk = fopen("time_taken.d", "w");
   	fprintf(tk, "Elapsed time is %lf seconds\n", (double)(toc-tic)/CLOCKS_PER_SEC);
   	fclose(tk);
	//printf("Value of the stoc constant is %le", -sqrt(2*alpha*KB_TIMES_TEMP)/sqrt(gama*ms*vol));
	
	
	/*if(shutdown=='y' || shutdown=='Y')
	system("shutdown -h +1");
	
	else;*/
	return;
	
}



void displayPrompt()
{
	char *s1 = "1. Heff only";
	char *s2 = "2. Heff and damping";
	char *s3 = "3. Heff, damping and STT";
	char *s4 = "4. Heff, damping, STT and FLT";
	char *s5 = "5. Heff, damping, fluctuation";
	char *s6 = "6. Heff, damping, STT and fluctuation";
	char *s7 = "7. Heff, damping, STT, FLT and fluctuation";
	char *s8 = "8. Only Thermal Fluctuations";
	printf("Enter your choice\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n", s1, s2, s3, s4, s5, s6, s7, s8);
}

void heffOnly()
{	
	const double dt = time_step;
	alpha=0.0;
	//Temporary variables for the rk-2
	initializeK();
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3]={0.0};
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	mv_t *system_cartesian = malloc(NO_TIME_STEPS*sizeof(mv_t));
	mang_t system_polar;

	system_cartesian[0] = (mv_t){sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
	system_polar = (mang_t){theta, phi};
	//Creates filename
	char **filename = generateFileNames(1);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}	
	int j=1;
	while(j<NO_TIME_STEPS)
	{
		heffCalc(system_cartesian[j-1]);		
		aj_var = aj(system_polar);
		bj_var = beta*aj_var;		
		for(i=0; i<3; i++)
			T[i] = heff[i];
		//printf("T is %lf %lf %lf\n", T[0], T[1], T[2]);

 		kp[0] = phiDot(system_polar, T , Tp) ;
         	kt[0] = thetaDot(system_polar, T, Tp) ;

         	
         	mpolar =(mang_t){system_polar.theta + dt*(kt[0]),  system_polar.phi + dt*(kp[0]) };
         	 

         	kp[1] = phiDot(mpolar, T, Tp);
         	kt[1] = thetaDot(mpolar, T, Tp);

         	system_polar = (mang_t) {system_polar.theta + dt/2.*(kt[0]  + kt[1]), system_polar.phi + dt/2.*(kp[0]  + kp[1])};
		
         	 	
		if(system_polar.theta<0)
		{
			system_polar.theta =-system_polar.theta;	
			system_polar.phi = M_PI + system_polar.phi;
		}
		if(system_polar.theta>M_PI)
		{
			system_polar.theta = TWO_PI - system_polar.theta;
			system_polar.phi = system_polar.phi + M_PI;
		}
		
		system_cartesian[j] = (mv_t) {sin(system_polar.theta)*cos(system_polar.phi), sin(system_polar.theta)*sin(system_polar.phi), cos(system_polar.theta)};
		

		ti = ti+dt;
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;

	}
	ti=0;
	for(j=0; j<NO_TIME_STEPS;j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, system_cartesian[j].mx, system_cartesian[j].my, system_cartesian[j].mz);
		ti+=dt;
	}
	free(system_cartesian);
	fclose(fp);
}

void heffDamp()
{
	const double dt = time_step;
	//Temporary variables for the rk-2
	initializeK();
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	mv_t *system_cartesian = malloc(NO_TIME_STEPS*sizeof(mv_t));
	mang_t system_polar;

	system_cartesian[0] = (mv_t){sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
	system_polar = (mang_t){theta, phi};
	//Creates filename
	char **filename = generateFileNames(2);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}	
	int j=1;
	while(j<NO_TIME_STEPS)
	{
		heffCalc(system_cartesian[j-1]);		
		aj_var = aj(system_polar);
		bj_var = beta*aj_var;		
		for(i=0; i<3; i++)
		{
			T[i] = heff[i] ;
			Tp[i] = alpha*heff[i] ;
		}
		//printf("T is %lf %lf %lf\n", T[0], T[1], T[2]);

 		kp[0] = phiDot(system_polar, T , Tp) ;
         	kt[0] = thetaDot(system_polar, T, Tp) ;

         	
         	mpolar =(mang_t){system_polar.theta + dt*(kt[0]),  system_polar.phi + dt*(kp[0]) };
         	 

         	kp[1] = phiDot(mpolar, T, Tp);
         	kt[1] = thetaDot(mpolar, T, Tp);

         	system_polar = (mang_t) {system_polar.theta + dt/2.*(kt[0]  + kt[1]), system_polar.phi + dt/2.*(kp[0]  + kp[1])};
		
         	 	
		if(system_polar.theta<0)
		{
			system_polar.theta =-system_polar.theta;	
			system_polar.phi = M_PI + system_polar.phi;
		}
		if(system_polar.theta>M_PI)
		{
			system_polar.theta = TWO_PI - system_polar.theta;
			system_polar.phi = system_polar.phi + M_PI;
		}
		
		system_cartesian[j] = (mv_t) {sin(system_polar.theta)*cos(system_polar.phi), sin(system_polar.theta)*sin(system_polar.phi), cos(system_polar.theta)};
		

		ti = ti+dt;
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;

	}
	ti=0;
	for(j=0; j<NO_TIME_STEPS;j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, system_cartesian[j].mx, system_cartesian[j].my, system_cartesian[j].mz);
		ti+=dt;
	}
	free(system_cartesian);
	fclose(fp);
}

void heffDampStt()
{
	const double dt = time_step;
	//Temporary variables for the rk-2
	initializeK();
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	mv_t *system_cartesian = malloc(NO_TIME_STEPS*sizeof(mv_t));
	mang_t system_polar;

	system_cartesian[0] = (mv_t){sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
	system_polar = (mang_t){theta, phi};
	//Creates filename
	char **filename = generateFileNames(3);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}	
	int j=1;
	while(j<NO_TIME_STEPS)
	{
		heffCalc(system_cartesian[j-1]);		
		aj_var = aj(system_polar);
		bj_var = beta*aj_var;		
		for(i=0; i<3; i++)
		{
			T[i] = heff[i]  - alpha*aj_var*mf[i];
			Tp[i] = alpha*heff[i] + aj_var*mf[i];
		}
		//printf("T is %lf %lf %lf\n", T[0], T[1], T[2]);

 		kp[0] = phiDot(system_polar, T , Tp) ;
         	kt[0] = thetaDot(system_polar, T, Tp) ;

         	
         	mpolar =(mang_t){system_polar.theta + dt*(kt[0]),  system_polar.phi + dt*(kp[0]) };
         	 

         	kp[1] = phiDot(mpolar, T, Tp);
         	kt[1] = thetaDot(mpolar, T, Tp);

         	system_polar = (mang_t) {system_polar.theta + dt/2.*(kt[0]  + kt[1]), system_polar.phi + dt/2.*(kp[0]  + kp[1])};
		
         	 	
		if(system_polar.theta<0)
		{
			system_polar.theta =-system_polar.theta;	
			system_polar.phi = M_PI + system_polar.phi;
		}
		if(system_polar.theta>M_PI)
		{
			system_polar.theta = TWO_PI - system_polar.theta;
			system_polar.phi = system_polar.phi + M_PI;
		}
		
		system_cartesian[j] = (mv_t) {sin(system_polar.theta)*cos(system_polar.phi), sin(system_polar.theta)*sin(system_polar.phi), cos(system_polar.theta)};
		

		ti = ti+dt;
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;

	}
	ti=0;
	for(j=0; j<NO_TIME_STEPS;j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, system_cartesian[j].mx, system_cartesian[j].my, system_cartesian[j].mz);
		ti+=dt;
	}
	free(system_cartesian);
	fclose(fp);
}

void heffDampSttFlt()
{
	const double dt = time_step;
	//Temporary variables for the rk-2
	initializeK();
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	mv_t *system_cartesian = malloc(NO_TIME_STEPS*sizeof(mv_t));
	mang_t system_polar;

	system_cartesian[0] = (mv_t){sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
	system_polar = (mang_t){theta, phi};
	//Creates filename
	char **filename = generateFileNames(4);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}	
	int j=1;
	while(j<NO_TIME_STEPS)
	{
		heffCalc(system_cartesian[j-1]);		
		aj_var = aj(system_polar);
		bj_var = beta*aj_var;		
		for(i=0; i<3; i++)
		{
			T[i] = heff[i] + (bj_var - alpha*aj_var)*mf[i];
			Tp[i] = alpha*heff[i] + (aj_var+ alpha*bj_var)*mf[i];
		}
		//printf("T is %lf %lf %lf\n", T[0], T[1], T[2]);

 		kp[0] = phiDot(system_polar, T , Tp) ;
         	kt[0] = thetaDot(system_polar, T, Tp) ;

         	
         	mpolar =(mang_t){system_polar.theta + dt*(kt[0]),  system_polar.phi + dt*(kp[0]) };
         	 

         	kp[1] = phiDot(mpolar, T, Tp);
         	kt[1] = thetaDot(mpolar, T, Tp);

         	system_polar = (mang_t) {system_polar.theta + dt/2.*(kt[0]  + kt[1]), system_polar.phi + dt/2.*(kp[0]  + kp[1])};
		
         	 	
		if(system_polar.theta<0)
		{
			system_polar.theta =-system_polar.theta;	
			system_polar.phi = M_PI + system_polar.phi;
		}
		if(system_polar.theta>M_PI)
		{
			system_polar.theta = TWO_PI - system_polar.theta;
			system_polar.phi = system_polar.phi + M_PI;
		}
		
		system_cartesian[j] = (mv_t) {sin(system_polar.theta)*cos(system_polar.phi), sin(system_polar.theta)*sin(system_polar.phi), cos(system_polar.theta)};
		

		ti = ti+dt;
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;

	}
	ti=0;
	for(j=0; j<NO_TIME_STEPS;j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, system_cartesian[j].mx, system_cartesian[j].my, system_cartesian[j].mz);
		ti+=dt;
	}
	free(system_cartesian);
	fclose(fp);
}

void heffDampFluc()
{
	#warning If broken, change constant data type to register data type
	const double dt = time_step;
	initializeK();	
	//Temporary variables for the rk-2
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double skt[2];
	//k- term for stochastic phi
	double skp[2];
	//Temporary variable
	int i, j, k, iter;
	
	// declaring window size
	double windowAvg = 0, window_old_avg=0;
	int window_size  =  NO_TIME_STEPS/100;
	//Creating Movingwindow
	double **window = (double *)malloc(NO_ENSEMBLES * sizeof(double *));
	for(i =  0 ; i< NO_ENSEMBLES ; i++){
		window[i] = (double *) malloc (sizeof(double) * window_size );
	}
	double *windowAvg_for10 = (double *) malloc (sizeof(double) * NO_TIME_STEPS);
	double *mz_for_10 = (double *) malloc (sizeof(double) * NO_TIME_STEPS);
	
	//Following are the variables for the ensemble calculation. variable names are self explainatory
	mv_t *ensemble_cartesian = malloc(NO_ENSEMBLES*sizeof(mv_t));

	mang_t *ensemble_polar = malloc(NO_ENSEMBLES*sizeof(mang_t));

	int *inversion_flag=malloc(NO_ENSEMBLES*sizeof(int));

	double *inversion_time = malloc(NO_ENSEMBLES*sizeof(double));

	mv_t *ensemble_average = malloc(NO_TIME_STEPS*sizeof(mv_t));
	//Creates filename
	char **filename = generateFileNames(5);

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;
	}	
	// file pointer for window average
	FILE *fp_window = fopen("windowFor10.d", "w");
	double ct = cos(theta), st = sin(theta);
	double cp = cos(phi), sp = sin(phi);
	
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		ensemble_polar[k] = ((mang_t){theta, phi});
		ensemble_cartesian[k] = ((mv_t){st*cp, st*sp, ct});
		inversion_time[k] = dt*NO_TIME_STEPS;
	}
	ensemble_average[0] = ensemble_cartesian[0];
	
	for(i = 0 ; i<NO_ENSEMBLES ; i++){
		for(j = 0 ; j<window_size ; j++){
			window[i][j] = ct;
		}
	}
	
	printf("Function started\n");
	for(j=1; j<NO_TIME_STEPS; j++)
	{
		
		for(k=0; k<NO_ENSEMBLES; k++)
		{
		
			heffCalc(ensemble_cartesian[k]);			
			aj_var = aj(ensemble_polar[k]);
			bj_var = beta*aj_var;		
			for(i=0; i<3; i++)
			{
				T[i] = heff[i];
				Tp[i] = alpha*heff[i];
			}//loop i ends

 			kp[0] = phiDot(ensemble_polar[k], T , Tp) ;
         		kt[0] = thetaDot(ensemble_polar[k], T, Tp) ;

			skp[0] = stocPhiDot(ensemble_polar[k], &seed);
         		skt[0] = stocThetaDot(ensemble_polar[k], &seed);

         	
         		mpolar = (mang_t){ensemble_polar[k].theta+dt*(kt[0])+ sqrt(dt)*skt[0], ensemble_polar[k].phi+dt*(kp[0])+sqrt(dt)*skp[0]};

         		kp[1] = phiDot(mpolar, T, Tp);
         		kt[1] = thetaDot(mpolar, T, Tp);
			

 			skp[1] = stocPhiDot(mpolar, &seed);
         		skt[1] = stocThetaDot(mpolar, &seed);
			//printf("Deterministic term: %e\n Stochastic term: %e\n", dt/2.*(kp[0]  + kp[1]), 1/2.*(sp[0] + sp[1])*sqrt(4*M_PI*ms*gama*dt));

         		ensemble_polar[k]  = (mang_t){ensemble_polar[k].theta + dt/2.*(kt[0]  + kt[1]) + 1/2.*(skt[0] + skt[1])*sqrt(dt), ensemble_polar[k].phi + dt/2.*(kp[0]  + kp[1]) + 1/2.*(skp[0] + skp[1])*sqrt(dt)};

			if(ensemble_polar[k].theta<0)
			{
				ensemble_polar[k].theta =-ensemble_polar[k].theta;	
				ensemble_polar[k].phi = M_PI + ensemble_polar[k].phi;
			}
			if(ensemble_polar[k].theta>M_PI)
			{
				ensemble_polar[k].theta =TWO_PI -ensemble_polar[k].theta;	
				ensemble_polar[k].phi = M_PI + ensemble_polar[k].phi;
			}
			if(ensemble_polar[k].phi>TWO_PI)
			{
				ensemble_polar[k].phi = ensemble_polar[k].phi - TWO_PI;
			}
			
			if(ensemble_polar[k].phi<0)
			{
				ensemble_polar[k].phi = ensemble_polar[k].phi + TWO_PI;
			}			
			
			ct = cos(ensemble_polar[k].theta); st = sin(ensemble_polar[k].theta);
			cp = cos(ensemble_polar[k].phi); sp = sin(ensemble_polar[k].phi);	
			if(j>=window_size)
			{
			        window_old_avg = 0;
			        for(iter = 0 ; iter < window_size ; iter++){
				        window_old_avg += window[k][iter];
			        }
			        window_old_avg /= window_size;

	        	}		
			ensemble_cartesian[k] = ((mv_t){st*cp, st*sp, ct});

			

			window[k][(j-1)%window_size] = ensemble_cartesian[k].mz;
			if(k==10)
			{
			        mz_for_10[j] = ensemble_cartesian[k].mz;
			}
			if(j>=window_size)
			{
			        windowAvg = 0;
			        for(iter = 0 ; iter < window_size ; iter++){
				        windowAvg += window[k][iter];
			        }
			        windowAvg /= window_size;
        			if(k == 10)
	        		{
	        			windowAvg_for10[j - window_size/2] = windowAvg;
	        			fprintf(fp_window,"%lf\t%lf\t%lf\t%e\t%lf\n",ti, windowAvg_for10[j - window_size/2], window_old_avg,fabs(windowAvg_for10[j - window_size/2]- window_old_avg), mz_for_10[j - window_size/2]);
	        		}
	        	}
			if(inversion_flag[k]==0 && fabs(windowAvg)>mz_tol )
			{
                                inversion_time[k]=ti;
                                inversion_flag[k]=1;
			}	        	
			

		
		}//loop k ends

		ensemble_average[j] = (mv_t) {0.0, 0.0, 0.0};
		

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			ensemble_average[j] = (mv_t){ensemble_average[j].mx + ensemble_cartesian[k].mx,ensemble_average[j].my + ensemble_cartesian[k].my , ensemble_average[j].mz + ensemble_cartesian[k].mz };
		
 		}
 		
 		ensemble_average[j] = (mv_t){ensemble_average[j].mx/(double)NO_ENSEMBLES, ensemble_average[j].my/(double)NO_ENSEMBLES, ensemble_average[j].mz/(double)NO_ENSEMBLES};

		
		ti = ti+dt;


		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));

	}//loop j ends
	ti=0.0;

	for(j=0; j<NO_TIME_STEPS; j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, ensemble_average[j].mx, ensemble_average[j].my, ensemble_average[j].mz);
		ti = ti+dt;
	}

	storeInversionTime(inversion_time, *(filename+1));
	free(ensemble_polar);
	free(ensemble_cartesian);
	free(ensemble_average);
	free(inversion_flag);
	free(inversion_time);
	fclose(fp);
	
	//closing window pointer 
	fclose(fp_window);

}
void heffDampSttFluc()
{
	#warning If broken, change constant data type to register data type
	const double dt = time_step;
	initializeK();	
	//Temporary variables for the rk-2
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double skt[2];
	//k- term for stochastic phi
	double skp[2];
	//Temporary variable
	int i, j, k, iter;
	
	// declaring window size
	double windowAvg = 0, window_old_avg=0;
	int window_size  =  NO_TIME_STEPS/100;
	//Creating Movingwindow
	double **window = (double *)malloc(NO_ENSEMBLES * sizeof(double *));
	for(i =  0 ; i< NO_ENSEMBLES ; i++){
		window[i] = (double *) malloc (sizeof(double) * window_size );
	}
	double *windowAvg_for10 = (double *) malloc (sizeof(double) * NO_TIME_STEPS);
	double *mz_for_10 = (double *) malloc (sizeof(double) * NO_TIME_STEPS);
	
	//Following are the variables for the ensemble calculation. variable names are self explainatory
	mv_t *ensemble_cartesian = malloc(NO_ENSEMBLES*sizeof(mv_t));

	mang_t *ensemble_polar = malloc(NO_ENSEMBLES*sizeof(mang_t));

	int *inversion_flag=malloc(NO_ENSEMBLES*sizeof(int));

	double *inversion_time = malloc(NO_ENSEMBLES*sizeof(double));

	mv_t *ensemble_average = malloc(NO_TIME_STEPS*sizeof(mv_t));
	//Creates filename
	char **filename = generateFileNames(6);

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;
	}	
	// file pointer for window average
	FILE *fp_window = fopen("windowFor10.d", "w");
	double ct = cos(theta), st = sin(theta);
	double cp = cos(phi), sp = sin(phi);
	
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		ensemble_polar[k] = ((mang_t){theta, phi});
		ensemble_cartesian[k] = ((mv_t){st*cp, st*sp, ct});
		inversion_time[k] = dt*NO_TIME_STEPS;
	}
	ensemble_average[0] = ensemble_cartesian[0];
	
	for(i = 0 ; i<NO_ENSEMBLES ; i++){
		for(j = 0 ; j<window_size ; j++){
			window[i][j] = ct;
		}
	}
	
	printf("Function started\n");
	for(j=1; j<NO_TIME_STEPS; j++)
	{
		
		for(k=0; k<NO_ENSEMBLES; k++)
		{
		
			heffCalc(ensemble_cartesian[k]);			
			aj_var = aj(ensemble_polar[k]);
			bj_var = beta*aj_var;		
			for(i=0; i<3; i++)
			{
				T[i] = heff[i]  - alpha*aj_var*mf[i];
				Tp[i] = alpha*heff[i] + aj_var*mf[i];
			}//loop i ends

 			kp[0] = phiDot(ensemble_polar[k], T , Tp) ;
         		kt[0] = thetaDot(ensemble_polar[k], T, Tp) ;

			skp[0] = stocPhiDot(ensemble_polar[k], &seed);
         		skt[0] = stocThetaDot(ensemble_polar[k], &seed);

         	
         		mpolar = (mang_t){ensemble_polar[k].theta+dt*(kt[0])+ sqrt(dt)*skt[0], ensemble_polar[k].phi+dt*(kp[0])+sqrt(dt)*skp[0]};

         		kp[1] = phiDot(mpolar, T, Tp);
         		kt[1] = thetaDot(mpolar, T, Tp);
			

 			skp[1] = stocPhiDot(mpolar, &seed);
         		skt[1] = stocThetaDot(mpolar, &seed);
			//printf("Deterministic term: %e\n Stochastic term: %e\n", dt/2.*(kp[0]  + kp[1]), 1/2.*(sp[0] + sp[1])*sqrt(4*M_PI*ms*gama*dt));

         		ensemble_polar[k]  = (mang_t){ensemble_polar[k].theta + dt/2.*(kt[0]  + kt[1]) + 1/2.*(skt[0] + skt[1])*sqrt(dt), ensemble_polar[k].phi + dt/2.*(kp[0]  + kp[1]) + 1/2.*(skp[0] + skp[1])*sqrt(dt)};

			if(ensemble_polar[k].theta<0)
			{
				ensemble_polar[k].theta =-ensemble_polar[k].theta;	
				ensemble_polar[k].phi = M_PI + ensemble_polar[k].phi;
			}
			if(ensemble_polar[k].theta>M_PI)
			{
				ensemble_polar[k].theta =TWO_PI -ensemble_polar[k].theta;	
				ensemble_polar[k].phi = M_PI + ensemble_polar[k].phi;
			}
			if(ensemble_polar[k].phi>TWO_PI)
			{
				ensemble_polar[k].phi = ensemble_polar[k].phi - TWO_PI;
			}
			
			if(ensemble_polar[k].phi<0)
			{
				ensemble_polar[k].phi = ensemble_polar[k].phi + TWO_PI;
			}			
			
			ct = cos(ensemble_polar[k].theta); st = sin(ensemble_polar[k].theta);
			cp = cos(ensemble_polar[k].phi); sp = sin(ensemble_polar[k].phi);	
			if(j>=window_size)
			{
			        window_old_avg = 0;
			        for(iter = 0 ; iter < window_size ; iter++){
				        window_old_avg += window[k][iter];
			        }
			        window_old_avg /= window_size;

	        	}		
			ensemble_cartesian[k] = ((mv_t){st*cp, st*sp, ct});

			

			window[k][(j-1)%window_size] = ensemble_cartesian[k].mz;
			if(k==10)
			{
			        mz_for_10[j] = ensemble_cartesian[k].mz;
			}
			if(j>=window_size)
			{
			        windowAvg = 0;
			        for(iter = 0 ; iter < window_size ; iter++){
				        windowAvg += window[k][iter];
			        }
			        windowAvg /= window_size;
        			if(k == 10)
	        		{
	        			windowAvg_for10[j - window_size/2] = windowAvg;
	        			fprintf(fp_window,"%lf\t%lf\t%lf\t%e\t%lf\n",ti, windowAvg_for10[j - window_size/2], window_old_avg,fabs(windowAvg_for10[j - window_size/2]- window_old_avg), mz_for_10[j - window_size/2]);
	        		}
	        	}
			if(inversion_flag[k]==0 && fabs(windowAvg)>mz_tol )
			{
                                inversion_time[k]=ti;
                                inversion_flag[k]=1;
			}	        	

		
		}//loop k ends

		ensemble_average[j] = (mv_t) {0.0, 0.0, 0.0};
		

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			ensemble_average[j] = (mv_t){ensemble_average[j].mx + ensemble_cartesian[k].mx,ensemble_average[j].my + ensemble_cartesian[k].my , ensemble_average[j].mz + ensemble_cartesian[k].mz };
		
 		}
 		
 		ensemble_average[j] = (mv_t){ensemble_average[j].mx/(double)NO_ENSEMBLES, ensemble_average[j].my/(double)NO_ENSEMBLES, ensemble_average[j].mz/(double)NO_ENSEMBLES};

		
		ti = ti+dt;


		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	

	}//loop j ends
	ti=0.0;

	for(j=0; j<NO_TIME_STEPS; j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, ensemble_average[j].mx, ensemble_average[j].my, ensemble_average[j].mz);
		ti = ti+dt;
	}

	storeInversionTime(inversion_time, *(filename+1));
	free(ensemble_polar);
	free(ensemble_cartesian);
	free(ensemble_average);
	free(inversion_flag);
	free(inversion_time);
	fclose(fp);
	
	//closing window pointer 
	fclose(fp_window);


}
void heffDampSttFltFluc()
{	
	#warning If broken, change constant data type to register data type
	const double dt = time_step;
	initializeK();	
	//Temporary variables for the rk-2
	mang_t mpolar;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double skt[2];
	//k- term for stochastic phi
	double skp[2];
	//Temporary variable
	int i, j, k, iter;
	
	// declaring window size
	double windowAvg = 0, window_old_avg=0;
	int window_size  =  NO_TIME_STEPS/500;
	//Creating Movingwindow
	double **window = (double *)malloc(NO_ENSEMBLES * sizeof(double *));
	for(i =  0 ; i< NO_ENSEMBLES ; i++){
		window[i] = (double *) malloc (sizeof(double) * window_size );
	}
	double *windowAvg_for10 = (double *) malloc (sizeof(double) * NO_TIME_STEPS);
	double *mz_for_10 = (double *) malloc (sizeof(double) * NO_TIME_STEPS);
	
	//Following are the variables for the ensemble calculation. variable names are self explainatory
	mv_t *ensemble_cartesian = malloc(NO_ENSEMBLES*sizeof(mv_t));

	mang_t *ensemble_polar = malloc(NO_ENSEMBLES*sizeof(mang_t));

	int *inversion_flag=malloc(NO_ENSEMBLES*sizeof(int));

	double *inversion_time = malloc(NO_ENSEMBLES*sizeof(double));

	mv_t *ensemble_average = malloc(NO_TIME_STEPS*sizeof(mv_t));
	//Creates filename
	char **filename = generateFileNames(7);

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;
	}	
	// file pointer for window average
	FILE *fp_window = fopen("windowFor10.d", "w");
	double ct = cos(theta), st = sin(theta);
	double cp = cos(phi), sp = sin(phi);
	
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		ensemble_polar[k] = ((mang_t){theta, phi});
		ensemble_cartesian[k] = ((mv_t){st*cp, st*sp, ct});
		inversion_time[k] = dt*NO_TIME_STEPS;
	}
	ensemble_average[0] = ensemble_cartesian[0];
	
	for(i = 0 ; i<NO_ENSEMBLES ; i++){
		for(j = 0 ; j<window_size ; j++){
			window[i][j] = ct;
		}
	}
	
	printf("Function started\n");
	for(j=1; j<NO_TIME_STEPS; j++)
	{
		
		for(k=0; k<NO_ENSEMBLES; k++)
		{
		
			heffCalc(ensemble_cartesian[k]);			
			aj_var = aj(ensemble_polar[k]);
			bj_var = beta*aj_var;		
			for(i=0; i<3; i++)
			{
				T[i] = heff[i] + (bj_var - alpha*aj_var)*mf[i];
				Tp[i] = alpha*heff[i] + (aj_var+ alpha*bj_var)*mf[i];
			}//loop i ends
			//printf("Damping torque: %lf. Alignment torque :%lf. ratio = %lf\n", 4*M_PI*ms*alpha*heff[2], 4*M_PI*ms*(aj_var+ alpha*bj_var), alpha*heff[2]/(aj_var+ alpha*bj_var) );
			//printf("Equillibrium vector : (%lf %lf %lf)\n", Tp[0], Tp[1], Tp[2]);

 			kp[0] = phiDot(ensemble_polar[k], T , Tp) ;
         		kt[0] = thetaDot(ensemble_polar[k], T, Tp) ;

			skp[0] = stocPhiDot(ensemble_polar[k], &seed);
         		skt[0] = stocThetaDot(ensemble_polar[k], &seed);

         	
         		mpolar = (mang_t){ensemble_polar[k].theta+dt*(kt[0])+ sqrt(dt)*skt[0], ensemble_polar[k].phi+dt*(kp[0])+sqrt(dt)*skp[0]};

         		kp[1] = phiDot(mpolar, T, Tp);
         		kt[1] = thetaDot(mpolar, T, Tp);
			

 			skp[1] = stocPhiDot(mpolar, &seed);
         		skt[1] = stocThetaDot(mpolar, &seed);
			//printf("Deterministic term: %e\n Stochastic term: %e\n", dt/2.*(kp[0]  + kp[1]), 1/2.*(sp[0] + sp[1])*sqrt(4*M_PI*ms*gama*dt));

         		ensemble_polar[k]  = (mang_t){ensemble_polar[k].theta + dt/2.*(kt[0]  + kt[1]) + 1/2.*(skt[0] + skt[1])*sqrt(dt), ensemble_polar[k].phi + dt/2.*(kp[0]  + kp[1]) + 1/2.*(skp[0] + skp[1])*sqrt(dt)};

			if(ensemble_polar[k].theta<0)
			{
				ensemble_polar[k].theta =-ensemble_polar[k].theta;	
				ensemble_polar[k].phi = M_PI + ensemble_polar[k].phi;
			}
			if(ensemble_polar[k].theta>M_PI)
			{
				ensemble_polar[k].theta =TWO_PI -ensemble_polar[k].theta;	
				ensemble_polar[k].phi = M_PI + ensemble_polar[k].phi;
			}
			if(ensemble_polar[k].phi>TWO_PI)
			{
				ensemble_polar[k].phi = ensemble_polar[k].phi - TWO_PI;
			}
			
			if(ensemble_polar[k].phi<0)
			{
				ensemble_polar[k].phi = ensemble_polar[k].phi + TWO_PI;
			}			
			
			ct = cos(ensemble_polar[k].theta); st = sin(ensemble_polar[k].theta);
			cp = cos(ensemble_polar[k].phi); sp = sin(ensemble_polar[k].phi);	
			if(j>=window_size)
			{
			        window_old_avg = 0;
			        for(iter = 0 ; iter < window_size ; iter++){
				        window_old_avg += window[k][iter];
			        }
			        window_old_avg /= window_size;

	        	}		
			ensemble_cartesian[k] = ((mv_t){st*cp, st*sp, ct});

			

			window[k][(j-1)%window_size] = ensemble_cartesian[k].mz;
			if(k==10)
			{
			        mz_for_10[j] = ensemble_cartesian[k].mz;
			}
			if(j>=window_size)
			{
			        windowAvg = 0;
			        for(iter = 0 ; iter < window_size ; iter++){
				        windowAvg += window[k][iter];
			        }
			        windowAvg /= window_size;
        			if(k == 10)
	        		{
	        			windowAvg_for10[j - window_size/2] = windowAvg;
	        			fprintf(fp_window,"%lf\t%lf\t%lf\t%e\t%lf\n",ti, windowAvg_for10[j - window_size/2], window_old_avg,fabs(windowAvg_for10[j - window_size/2]- window_old_avg), mz_for_10[j - window_size/2]);
	        		}
	        	}
			
			/*****************************THIS IS WHERE YOU NEED TO PUT THE NEW SWITCHING TIME DEFINITION***************/
			if(inversion_flag[k]==0 && fabs(windowAvg)>mz_tol )
			{
                                inversion_time[k]=ti;
                                inversion_flag[k]=1;
			}

		
		}//loop k ends

		ensemble_average[j] = (mv_t) {0.0, 0.0, 0.0};
		

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			ensemble_average[j] = (mv_t){ensemble_average[j].mx + ensemble_cartesian[k].mx,ensemble_average[j].my + ensemble_cartesian[k].my , ensemble_average[j].mz + ensemble_cartesian[k].mz };
		
 		}
 		
 		ensemble_average[j] = (mv_t){ensemble_average[j].mx/(double)NO_ENSEMBLES, ensemble_average[j].my/(double)NO_ENSEMBLES, ensemble_average[j].mz/(double)NO_ENSEMBLES};

		
		ti = ti+dt;


		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));

	}//loop j ends
	ti=0.0;

	for(j=0; j<NO_TIME_STEPS; j++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", ti, ensemble_average[j].mx, ensemble_average[j].my, ensemble_average[j].mz);
		ti = ti+dt;
	}

	storeInversionTime(inversion_time, *(filename+1));
	free(ensemble_polar);
	free(ensemble_cartesian);
	free(ensemble_average);
	free(inversion_flag);
	free(inversion_time);
	fclose(fp);
	
	//closing window pointer 
	fclose(fp_window);
}

void flucOnly()
{
	/*register double dt = time_step;	
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double st[2];
	//k- term for stochastic phi
	double sp[2];
	//Temporary variable
	int i;
	
	double mx[NO_ENSEMBLES], my[NO_ENSEMBLES], mz[NO_ENSEMBLES];
	double mxavg, myavg, mzavg;
	//Creates filename
	char **filename = generateFileNames(8);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	int j, k;

	double theta_en[NO_ENSEMBLES], phi_en[NO_ENSEMBLES];
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		theta_en[k] = theta;
		phi_en[k] = phi;

	}

	for(j=1; j<NO_TIME_STEPS; j++)
	{
		
		for(k=0; k<NO_ENSEMBLES; k++)
		{

			


			sp[0] = stocPhiDot(theta_en[k], phi_en[k], &seed);
         		st[0] = stocThetaDot(theta_en[k], phi_en[k], &seed);

         	
         		mphi = phi_en[k] + sqrt(dt)*sp[0] ;
         		mtheta = theta_en[k] + sqrt(dt)*st[0] ;
			
			sp[1] = stocPhiDot(mtheta, mphi, &seed);
         		st[1] = stocThetaDot(mtheta, mphi, &seed);

         		phi_en[k] = phi_en[k] +  (1/2.)*(sp[0] + sp[1])*sqrt(dt);
         		theta_en[k] = theta_en[k] + (1/2.)*(st[0] + st[1])*sqrt(dt);

			if(theta_en[k]<0)
			{
				theta_en[k] =-theta_en[k];	
				phi_en[k] = M_PI + phi_en[k];
			}
			if(theta_en[k]>M_PI)
			{
				theta_en[k] = TWO_PI - theta_en[k];
				phi_en[k] = phi_en[k] + M_PI;
			}	
		
			mx[k] = sin(theta_en[k])*cos(phi_en[k]);
			my[k] = sin(theta_en[k])*sin(phi_en[k]);
			mz[k] = cos(theta_en[k]);
		}

		mxavg = 0;
 		myavg = 0;
 		mzavg = 0;

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			mxavg+=mx[k];
 			myavg+=my[k];
 			mzavg+=mz[k];
 		}
 		
 		mxavg = mxavg/NO_ENSEMBLES;
 		myavg = myavg/NO_ENSEMBLES;
 		mzavg = mzavg/NO_ENSEMBLES;
		

		ti = ti+dt;
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mxavg, myavg, mzavg);
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));

	}
	fclose(fp);*/

}

void initializeK()
{
	K = 1/(1+ alpha*alpha);
}


double thetaDot(mang_t polar_coord, double *T, double *Tp)
{

	static double ct, st, cp, sp, val;
	ct = cos(polar_coord.theta); st = sin(polar_coord.theta);
	cp = cos(polar_coord.phi); sp = sin(polar_coord.phi);
	val= ( -K*(( sp*T[0] ) - (cp*T[1]) - (cp*ct*Tp[0]) - (sp*ct*Tp[1]) + (st*Tp[2]) ));
	//printf("Theta dot is %lf\n", val);
	return val;
}




double phiDot(mang_t polar_coord, double *T, double *Tp)
{
	static double ct, st, cp, sp, val;
	ct = cos(polar_coord.theta); st = sin(polar_coord.theta);
	cp = cos(polar_coord.phi); sp = sin(polar_coord.phi);
	val = ( -K*( (cp*ct*T[0]/st) + (T[1]*ct*sp/st)  - T[2] + (Tp[0]*sp/st) - (Tp[1]*cp/st) ));	
	return val;
}

double stocThetaDot(mang_t polar_coord, long *seed)
{
	static double ct, st, cp, sp, val;
	ct = cos(polar_coord.theta); st = sin(polar_coord.theta);
	cp = cos(polar_coord.phi); sp = sin(polar_coord.phi);	
	val = gauss(seed)*(sp - alpha*ct*cp) - gauss(seed)*(cp - alpha*ct*sp) + alpha*gauss(seed)*st;
	val =  -K*sqrt(2*alpha*KB_TIMES_TEMP)/sqrt(gama*ms*vol)*val;
	//printf("Hlf is %.4f\n", val/4*M_PI*ms);
	return val/(4*M_PI*ms)*sqrt(4*M_PI*ms*gama);
	

}

double stocPhiDot(mang_t polar_coord, long *seed)
{
	static double ct, st, cp, sp, val;
	ct = cos(polar_coord.theta); st = sin(polar_coord.theta);
	cp = cos(polar_coord.phi); sp = sin(polar_coord.phi);	
	val = gauss(seed)*(ct*cp + alpha*sp) + gauss(seed)*(ct*sp - alpha*cp) - gauss(seed)  ;
	val = -K*sqrt(2*alpha*KB_TIMES_TEMP)/sqrt(gama*ms*vol)*val/st;
	val =  val/(4*M_PI*ms)*sqrt(4*M_PI*ms*gama);
	//printf("Value of stocPhiDot is %e\n", val);
	return val;
	
}
double aj(mang_t polar_coord)
{
	static double dot_product, g, value, ct, st, cp, sp;
	ct = cos(polar_coord.theta); st = sin(polar_coord.theta);
	cp = cos(polar_coord.phi); sp = sin(polar_coord.phi);
	dot_product  = st*cp*mf[0] + st*sp*mf[1] + ct*mf[2];	
	g = eta/(1 + 0.1*dot_product);
	value = (HCROSS*g*current)/(2*ECHARGE*4*M_PI*ms*vol);
	value = value/(4*M_PI*ms);
	return value;

}

void heffCalc(mv_t cartesian_coord)
{
	heff[0] = hext[0] - (hde[0] - han[0])*cartesian_coord.mx;
	heff[1] = hext[1] - (hde[1] - han[1])*cartesian_coord.my;
	heff[2] = hext[2] - (hde[2] - han[2])*cartesian_coord.mz;
}


void storeInversionTime(double *inversion_time, char *filename)
{
	
	FILE *fp = fopen(filename, "w");
	int i;
	for(i=0; i<NO_ENSEMBLES; i++)
	{
		fprintf(fp, "%d\t%lf\n", i+1, inversion_time[i]);

	}
	fclose(fp);
}





