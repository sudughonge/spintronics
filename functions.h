/*
 * Simulation for Spintronics and Spin Transfer Torque
 * Code writen by Sudarshan Ghonge and Ajinkya Gavane
 * Department of Physics,
 * Birla Institute Institute of Technology and Science, Pilani. Hyderabad Campus
 *
 * 
 * The following is a program to simulate the dynamics of a soft iron layer in Magnetic Tunneling Junction.
 * The Landu-Lifhittz-Gilbert Equation has been used.
 * Refer http://en.wikipedia.org/wiki/Landau%E2%80%93Lifshitz%E2%80%93Gilbert_equation
 
 * The terms used are 
 * 1. Effective Magnetic Field
 * 2. Damping
 * 3. Spin Transfer Torque
 * 4. Field Like Toruqe
 * 5. Stochastic Field simulated using white Gaussian Noise

 * The load of solving for three quantities namely mx, my and mz -the three components of the magnetization of the soft iron magnet
 * has been reduced by assuming that none of the terms lead to dissipation and the magnitude remains constant and hence the dynamics of only the
 * angular variables has been computed.

 * The dynamics have been done on the time variable using Runge-Kutta 2 method.
 

 */


#include<stdio.h>

#include<time.h>
#include<string.h>
#define HCROSS 1.054e-27   //The reduced planck's constant in cgs units
#define ECHARGE 1.602e-20  //The electric charge in emu
#define KB  1.380e-16 //Value of "Kb", the Boltzmann Constant in cgs
#define TWO_PI 2*M_PI   //Value of 2*pi
#define CLEARBUF() char ch; while(ch = getchar() != '\n' && ch!= EOF); //CLEARBUF is used to clear the input buffer. 
#define TIME_STEP 1e-13

/*Structure of the magnetic moment of the soft ferromagnet-the mvector in the cartesian co-ordinates.
Typedef as 'mv_t'*/
typedef struct MVector
{
	double mx;
	double my;
	double mz;
}mv_t;

/* Structure of the angular co-ordinates of the soft ferromagnet-the mvector in spherical polar co-ordinates.
Typedef as 'mang_t'*/
typedef struct MAngle
{
	double theta;
	double phi;

}mang_t;

/*

	Declaration of static variables which retain values of the parameters and the function dostuff() makes them dimensionless

*/
extern double alpha, gama, current, ms, eta, beta, vol, K;
extern double mz_tol, mz_der_tol;
extern double heff[3], mf[3], hext[3], han[3], hde[3];
extern double theta, phi;
extern int NO_TIME_STEPS, NO_ENSEMBLES;
extern double KB_TIMES_TEMP;

extern double ti, time_step;

extern double aj_var, bj_var;


extern long seed;

/*
	Returns value of aj which is instantaneously dependent on m-hat
*/
double aj(mang_t);

void dostuff( double , double ,double , double , double , double , double , double *, double *, double *,  double *, double , double , int , double , int , double , double, double , double  );

/*
	Returns value of heff using current value of mx, my and mz
*/
void heffCalc(mv_t);

/*
	The following two functions return the dimensionless value of the deterministic part of the dynamics of theta and phi.
	The value returned depends on the passed value of theta, phi, T and Tp.
*/

double thetaDot(mang_t , double *, double *);
double phiDot(mang_t, double *, double *);
/*
	The following two functions return the dimensionless value of the stochastic part of the dynamics of theta and phi.
	The value returned depends on the passed value of theta, phi and on the random number generator.
*/

double stocThetaDot(mang_t , long *);
double stocPhiDot(mang_t, long *);

/*
	The following stores the inversion time and dwelling time of the soft magnet
*/

void storeInversionTime(double *, char*);
void storeDwellingTime(double *, char *);

/*
Initialzing K
*/
void initializeK();

char **generateFileNames(int );


/*
	Below are various modes of the simulation 
*/

/*
 * 1. Simulates the presence of only the Effective field
 * 2. Simulates the presence of the Effective field in the presence of Damping
 * 3. Simulates the presence of the Effective field, damping and Spin Transfer Torque
 * 4. Simulates the presence of all deterministic terms: Effective field, damping, Spin Transfer Torque and Field Like Torque
 * 5. Simulates the presence of the Effective field, damping and Fluctuating field
 * 6. Simulates the presence of the Effective field, damping, Spin Transfer Torque and Fluctuating field
 * 7. Simulates the presence of the Effective field, damping, STT, FLT and Fluctuating field
 * 8. Only the Fluctuating field.
 */

//1.
void heffOnly();
//2.
void heffDamp();
//3.
void heffDampStt();
//4.
void heffDampSttFlt();
//5
void heffDampFluc();
//6
void heffDampSttFluc();
//7
void heffDampSttFltFluc();
//8
void flucOnly();

void displayPrompt();


