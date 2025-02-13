/******************************************************************************
* C implementation of:                                                        *
* A Quantitative Model of Human Jejunal Smooth Muscle Cell Electrophysiology  *
*                                                                             *
* This draft version of the code for public release was prepared by:          *
* YC Poh (email address: poh.yongcheng@gmail.com)                             *
*                                                                             *
* Paper details:                                                              *
* Authors: Poh YC, Corrias A, Cheng N and Buist ML                            *
* PLoS ONE 7(8): e42385. doi:10.1371/journal.pone.0042385                     *
*                                                                             *
* Feedback & comments are welcomed                                            *
*                                                                             *
* The CellML & Matlab versions of this model are available for download from  *
* the Computational Bioengineering Laboratory (NUS) website:                  *
* http://www.bioeng.nus.edu.sg/compbiolab/                                    *
******************************************************************************/

//List of all subroutines/code obtained from
//"Numerical Recipes in C - The Art of Scientific Computing (2nd Ed)"
//used to solve the hJSMC model:
//1. ludcmp.c
//2. lubksb.c
//3. nrutil.c
//4. nr.h
//5. nrutil.h
//
//NOTE that the aforementioned Numerical Recipes subroutines are not provided
//in this download package due to copyright considerations.
//
//To get the Numerical Recipes code into your computer, you can copy them from
//the Numerical Recipes book itself: http://apps.nrbook.com/c/index.html
//
//Do also refer to the online NR-book for further information, or
//feel free to contact us to get help with implementing this hJSMC model (in C). 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h" //via backward Euler & LU decomposition

//Prototype
void ludcmp(double **a, int n, int *indx, double *d); //Numerical Recipes
void lubksb(double **a, int n, int *indx, double b[]); //Numerical Recipes
void norm_stat_val(int n, double b[]);

int main()
{

/* -------------------- */
/* Variable declaration */
/* -------------------- */

//Use enum for the states
enum Na_states
{
    ONa=1,C1=2,C2=3,C3=4,I1=5,I2=6
};

enum CaL_states
{
  C0_CaL=1,C1_CaL=2,C2_CaL=3,C3_CaL=4,O_CaL=5,Ivf_CaL=6,Ivs_CaL=7,C0Ca_CaL=8,C1Ca_CaL=9,C2Ca_CaL=10,C3Ca_CaL=11,ICa_CaL=12,IvfCa_CaL=13,IvsCa_CaL=14
};

enum BK_states
{
  O0_BK=1,O1_BK=2,O2_BK=3,O3_BK=4,O4_BK=5,C0_BK=6,C1_BK=7,C2_BK=8,C3_BK=9,C4_BK=10
};

//Define cell variables and constants
double t,V,dVdt,R,T,F,Cm,Am,Vcell;//time;membrane voltage;change in membrane voltage over time;ideal gas constant;model temperature;Faraday constant;membrane capacitance;volume of cell

//Define control & misc. variables
double Texp; //experimental temperature
double dt; //time step size
int printfidx; //index for printing results into file
int printsidx; //index for printing results onto screen
int printfrequency; //the frequency of printing results
int i; //index for selective printing, runs parallel with time steps

//Define variables for L-type calcium current
double a,a0,a1,a2,a3,b,b0,b1,b2,b3,yf,ys,pf,ps,lf,ls,wf,ws,wsf,wfs,delta,theta; //rate variables
double ICaL, GCaL; //CaL current and conductance
double CalciumDependenceOn; //1 is to switch off calcium dependency; 0 to switch off
double *B_CaL; //NR B vector
double **A_CaL; //NR A matrix
double d_nr_CaL; //NR
int *indx_CaL; //NR
int i_nr_CaL,n_nr_CaL; //NR


//Define variables for T-type calcium current
double ICaT, GCaT; //CaT current and conductance
double dCaT, fCaT; //inactivation and activation gating variables
double dCaT_inf, fCaT_inf; //steady-state gating variables
double tau_dCaT, tau_fCaT; //time constant taus of gating variables
double ddCaTdt, dfCaTdt; //time derivative of gating variables

//Define variables for voltage-dependent potassium current
double IKv, GKv; //Kv current and conducance
double xKv, yKv; //inactivation and activation gating variables
double xKv_inf, yKv_inf; //steady-state gating variables
double tau_xKv, tau_yKv; //time constant taus of gating variables
double dxKvdt, dyKvdt; //time derivative of gating variables

//Define variables for voltage-dependent and calcium-sensitive potassium current
double IBK, GBK; //BK current and conductance
double qCO,commonval1,k_C0O0,k_C1O1,k_C2O2,k_C3O3,k_C4O4; //rate variables
double qOC,commonval2,k_O0C0,k_O1C1,k_O2C2,k_O3C3,k_O4C4; //rate variables
double k_ON,k_OFF_C,k_OFF_O; //rate variables
double k_C0C1,k_C1C2,k_C2C3,k_C3C4,k_C4C3,k_C3C2,k_C2C1,k_C1C0; //rate variables
double k_O0O1,k_O1O2,k_O2O3,k_O3O4,k_O4O3,k_O3O2,k_O2O1,k_O1O0; //rate variables
double F_BK,R_BK,T_BK,V_m,Ca_free,time_convert; //Faraday's constant; ideal gas constant;temperature,membrane voltage;free intracellular calcium;time conversion --> for use in the BK model
double paramBK[13]={21.6200042030946,0.868781002082468,0.028063354826539,0.781455969444612,44.323856772419,318108.401217897,144173.645394299,32659.4033438573,95.3119836990264,0.105563199075352,40.6326276428765,1,1};
double *B_BK; //NR B vector
double **A_BK; //NR A matrix
double d_nr_BK; //NR
int *indx_BK; //NR
int i_nr_BK,n_nr_BK; //NR

//Define variables for sodium current
double k_O_I1,k_I1_I2,k_C3_C2,k_C2_C1,k_C1_O,k_I2_I1,k_C2_C3,k_C1_C2,k_O_C1,k_I1_C1,k_C1_I1,k_I1_O; //rate variables
double INa, GNa; //Na current and conductance
double parval[36]={0.307629865993612,0.00605347675235982,0.0514901233349612,-0.0468650997035773,-0.0691017278134239,0.00319450166548765,
-0.156598571141704,0.0583527292113407,0.00931925729714207,0.0410752869450889,2.67926519075303,0.00614677585983095,
-0.0990741811465619,0.0364411841506491,0.363520774026992,0.0771931771113635,-13.3349821299039,-0.252888485042226,
-2.48395798104711,0.0204061599230419,-0.0634382282818863,0.00466834101392052,1.61639957785747,0.0277351044118711,
0.000525484598535321,1.4496348286393,1.53287364507221,0.00392388356048677,0.554315057682815,3.15656890165025,
2.39152928164775,1.90461121472293,0.00021688112860783,0.120515381956888,-9.60280306718205,0.0830250234337321 };
double *B; //NR B vector
double **A;//NR A matrix
double d_nr; //NR
int *indx; //NR
int i_nr,n_nr; //NR

//Define variables for non-selective cationic current
double INSCC,GNSCC; //NSCC current and conductance
double mNSCC,rAch,hCa; //gating variable, aceytlcholine dependency, and calcium dependency of the NSCC channels
double mNSCC_inf,tau_mNSCC,dmNSCCdt; //steady-state value of gating variable, time constant tau of gating variable, time derivative of mNSCC gating variable

//Define variables for sodium potassium pump
double INaK,PNaK; //NaK current and maximum PNaK current
double numINaK,denomINaK; //numerator and denominator variables of INaK
double KmK,KmNa1; //Ko half-saturation constant of INaK and Nai half-saturation constant of INaK respectively

//Define variables for sodium calcium exchanger
double INCX,KNCX; //NCX current and maximum NCX current
double numINCX,denomINCX; //numerator and denominator variables of INCX
double KmNa2,KmCa,Ksat,gamma; // Nai half-saturation constant and Cai half-saturation constant respectively, and saturation factor and voltage dependence parameter respectively

//Define variables for the ICC stimlus current
double Num_slow_waves,Slow_wave_period,Init_equilibration;
double Complete_slow_wave_marker,ICC_peak_duration,ICC_plateau_duration;
double localt,kSlope;
double VICC,VICCrest,VICCamp,part1,part2,Gcouple,Istim;
double ICC_ftr1, ICC_ftr2;

//Define variables for concentration tracking
double dNaidt,dCaidt,dKidt; //time derivatives of ionic concentrations
double Nai,Nao,Cai,Cao,Ki,Ko; //ionic concentrations in and out of cell
double ENa,ECa,EK,ENSCC; //Nernst potential
double z; //calcium valence

//Define variables for temperature correction and Q10
double T_correction_Na, Q10Na; //For Na
double T_correction_Ca, Q10Ca; //For Ca
double T_correction_K, Q10K; //For K
double T_correction_CaL; //For CaL

//Define variables for calcium buffering
//"fitted" to reproduce physiological Cai values
double theta_cabuff, theta_1_cabuff, theta_2_cabuff;
double hill_coeff_CRT=1.0;
double hill_coeff_CaM=4.0;
double CRT_max=0.034;//mM
double CaM_max=0.012;//mM
double K_CRT=0.0009;//mM
double K_CaM=0.0001;//mM
double Cai_free;//mM

//Define variable for sensitivity analysis
double increase_var = 1.3;
double decrease_var = 0.7;

//NS for Na, K
double INS_Na, INS_K;
double gNSNa=0.022488;
double gNSK=0.017512;

/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* ------------ */
/* Define file  */
/* ------------ */

FILE *pResults= fopen ("test.xls","w"); // Declarations of files where to store results
fprintf( pResults, "time\tVm\tICaL\tICaT\tIKv\tIBK\tINa\tINCX\tINaK\tIstim\tVICC\tCai\tCai_free\tNai\tKi\n" ); //Print results into file

/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* ------- */
/* Control */
/* ------- */

CalciumDependenceOn=1;//For L-type calcium channels: 1, switch on Ca dependency; 0 switch it off
i=0;//index for selective printing, runs parallel with time steps
Num_slow_waves=180.0; //number of ICC slow waves
printfrequency=100*10; //frequency of results printing
dt=0.1; //msec

/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* ----------------------- */
/* Assign parameter values */
/* ----------------------- */

//Cell and general parameter values
T=310.0; //model temperature, K
Texp=297.0; //experiment temperature, K
R=8.314; //ideal gas constant, J/(molK)
Am=0.0077; //cell membrane surface area,mm^2
F=96.48534; //Faraday's constant, C/mmol
z=2.0; //calcium valence
Vcell=0.0000035; //cell volume, mm^3
Cm=50.0; //cell capacitance, pF

//ICC stimulus
//Values for Hwang (7.5cpm); rest=-67, amp=47, Gcouple=4
//values for Lee (6cpm); rest=-58, amp=32, Gcouple=2.6
Init_equilibration=0.0; //to equilibriate the system before proper ICC stimulation, msec
Complete_slow_wave_marker=Init_equilibration;

//Lee et al 6cpm
Slow_wave_period=10000.0;//msec
VICCrest=-58.0;//mV
VICCamp=32.0; //mV
Gcouple=2.6;//4;//nS

localt=0.0; //msec
ICC_peak_duration=300.0; //msec
ICC_plateau_duration=9700.0;//msec
ICC_ftr1=12000.0;
ICC_ftr2=300.0;
kSlope=450.0*ICC_ftr1/9000.0; //msec

GCaL=1.5*0.96;//nS
GNa=25.1;//nS
GBK=80.0;//nS
GCaT=0.0425;//nS
GKv=1.5*0.6811;//nS
KNCX=1992.1865;//pA
PNaK=1.75*0.85*2.8*3.88886015;//pA

//Q10
Q10Na=2.45;
Q10Ca=2.1;
Q10K=3.1;

/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */

/* ---------- */
/* Initialize */
/* ---------- */

//Initialize the NR variables of the Markov models
//CaL model
n_nr_CaL=14;
indx_CaL=ivector(1,n_nr_CaL);
B_CaL=dvector(1,n_nr_CaL);
A_CaL=dmatrix(1,n_nr_CaL,1,n_nr_CaL);
//BK model
n_nr_BK=10;
indx_BK=ivector(1,n_nr_BK);
B_BK=dvector(1,n_nr_BK);
A_BK=dmatrix(1,n_nr_BK,1,n_nr_BK);
//Na model
n_nr=6;
indx=ivector(1,n_nr);
B=dvector(1,n_nr);
A=dmatrix(1,n_nr,1,n_nr);

//Membrane potential
V=-59.0; //mV

//Initial ionic concentrations
Ki=150;//mM
Ko=5.4;//mM
Nai=10.0;//mM
Nao=140.0;//mM
Cai=54.6*0.9e-4;//mM; this is total Cai
Cao=2.0;//mM
Cai_free=1.4*0.9e-4;//mM; this is free Cai

//This segment is uncommented when looking at the SMC Vm w/0 ICC stim
////To match voltage clamp conditions of Farrugia et al's study on whole hJSMC outward currents
//T=297.0;//K
//Ki=150.0;//mM
//Ko=5.9;//mM

GKv=GKv;

//Temperature correction coefficients
T_correction_Na=1.0*pow(Q10Na,(T-Texp)/10.0);
T_correction_Ca=pow(Q10Ca,(T-Texp)/10.0);
T_correction_K=pow(Q10K,(T-Texp)/10.0); //affects only IKv
T_correction_CaL=1.0*pow(Q10Ca,(T-310.0)/10.0);//-310K 'cuz model was constructed with that temperature

//Initial Nernst potentials
ECa=(R*T/(2.0*F))*log(Cao/Cai_free); //mV
EK=(R*T/F)*log(Ko/Ki); //mV
ENa=(R*T/F)*log(Nao/Nai); //mV

//Initial values for L-type calcium
B_CaL[C0_CaL]= 0.563168;
B_CaL[C1_CaL] = 0.167247;
B_CaL[C2_CaL] = 0.0186255;
B_CaL[C3_CaL] = 0.000921887;
B_CaL[O_CaL] = 1.71111e-05;
B_CaL[Ivf_CaL] = 1.85565e-05;
B_CaL[Ivs_CaL]  = 4.68845e-06;
B_CaL[C0Ca_CaL] = 0.0202722;
B_CaL[C1Ca_CaL] = 0.003;
B_CaL[C2Ca_CaL] = 0.00602034;
B_CaL[C3Ca_CaL] = 0.000670459;
B_CaL[ICa_CaL] = 3.3185e-05;
B_CaL[IvfCa_CaL] = 6.15944e-07;
B_CaL[IvsCa_CaL] = 0.0;

//Initial values for T-type calcium
dCaT= 0.0202;
fCaT= 1.0;

//Initial values voltage-dependent potassium current
xKv= 0.2;
yKv= 0.99;

//Initial values for voltage-dependent and calcium-sensitive potassium channels
B_BK[O0_BK]=0.1;
B_BK[O1_BK]=0.1;
B_BK[O2_BK]=0.1;
B_BK[O3_BK]=0.1;
B_BK[O4_BK]=0.1;
B_BK[C0_BK]=0.1;
B_BK[C1_BK]=0.1;
B_BK[C2_BK]=0.1;
B_BK[C3_BK]=0.1;
B_BK[C4_BK]=0.1;

//Initial values for sodium current
B[ONa] = 1.25534e-05;
B[C1] = 0.0127485 ;
B[C2] = 0.0343829;
B[C3] = 0.136166;
B[I1] = 0.000258498;
B[I2] = 0.0252312;


/* --------------------------------------------------------------- */
/* --------------------------------------------------------------- */


/* --------------------------------------- */
/* ---------start of time loop------------ */
/* --------------------------------------- */
for (t=0;t<((Num_slow_waves*Slow_wave_period)+Init_equilibration);t=t+dt)
{

        /* ----------------- */
        /* Numerical methods */
        /* ----------------- */
        //Markov models solved using backward Euler method with LU decomposition
        //All else solved using forward Euler method

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ------------ */
        /* ICC stimulus */
        /* ------------ */

        //Regular stimulus
        if (t<=Complete_slow_wave_marker+ICC_plateau_duration+ICC_peak_duration && t>=Complete_slow_wave_marker)
        {
                if(t<=Complete_slow_wave_marker+ICC_peak_duration && t>=Complete_slow_wave_marker)
                {
                        VICC=VICCrest+VICCamp*(localt/ICC_ftr2);
                }
                else
                {
                        part1=1.0+exp(-ICC_ftr1/(2.0*kSlope));
                        part2=1.0/(1.0+exp((localt-ICC_ftr2-(ICC_ftr1/2.0))/kSlope));
                        VICC=VICCrest+VICCamp*part1*part2;
                }
        localt=localt+dt;
        }
        else
        {
                VICC=VICCrest;
        }

        if (t>=Complete_slow_wave_marker+Slow_wave_period)
        {
                Complete_slow_wave_marker=Complete_slow_wave_marker+Slow_wave_period;
                localt=0.0;
        }


//    //ICC stimulus inhibited by 2APB action
//    if (i==((140*10000)+Init_equilibration)/dt) //simulate effect of 50microM of 2-APB effect on ICC-MY Ca cycling --> lee et al; 10000 is the original slow wave period
//    {
//        Slow_wave_period=12245.0; //msec about 4.9cpm
//        VICCrest=-56; //-62;//mV
//        VICCamp=0.329*33.5; //mV. amplitude was about 32.9% of original b4 blockage
//        localt=0.0; //msec
//        ICC_peak_duration=1.2*300.0; //msec. time to peak was increased by about 1.2 times (from 3.25 to 3.89)
//        ICC_plateau_duration=Slow_wave_period-1.2*300.0; //5000.0; //msec
//        ICC_ftr1=12000;
//        ICC_ftr2=300;
//        kSlope=450.0*ICC_ftr1/9000; //msec
//
//        Complete_slow_wave_marker=0.0;
//        localt=0.0;
////        told=t;
////        t=0;
//    }


//    //Total inhibition of ICC stimulation
//    if (i>=((0*Slow_wave_period)+Init_equilibration)/dt)
//    {
//       Gcouple=0.0;
//    }

        Istim=Gcouple*(VICC-V); //pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */



        /* ------------------------------- */
        /* Non-selective cationic currents */
        /* ------------------------------- */
        INS_Na=gNSNa*(V-ENa);
        INS_K=gNSK*(V-EK);


        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ---------------------- */
        /* L-type calcium current */
        /* ---------------------- */

        //Rate equations
        a=T_correction_CaL*0.7310*exp(V/(30.0));//0.7322
        a0=4.0*a;
        a1=3.0*a;
        a2=2.0*a;
        a3=a;

        b=T_correction_CaL*0.2149*exp(-V/(40.0));//0.2144
        b0=b;
        b1=2.0*b;
        b2=3.0*b;
        b3=4.0*b;

        yf=T_correction_CaL*0.4742*exp(V/(10.0));//0.4691
        ys=T_correction_CaL*0.005956*exp(-V/(40.0));//0.0059

        pf=T_correction_CaL*0.02197*exp(V/(500.0));//0.0217
        ps=T_correction_CaL*0.002320*exp(-V/(280.0));//0.0023

        lf=T_correction_CaL*0.01407*exp(-V/(300.0));//0.0141
        ls=T_correction_CaL*0.01213*exp(V/(500.0));//0.0121

        wf=(b3*lf*yf)/(a3*pf);
        ws=(b3*ls*ys)/(a3*ps);
        wsf=(ls*pf)/lf;
        wfs=ps;

        if (CalciumDependenceOn)
        {
                delta=T_correction_CaL*4.0/(1.0+1.0/(Cai_free));//calcium dependency in mM
                theta=T_correction_CaL*0.01;
        }
        else
        {
                delta=0.0;
                theta=0.0;
         }

        //Equation system rearranged into Ax=b format, following backward Euler implementation
        //Assign the row vectors of the A matrix

        A_CaL[C0_CaL][C0_CaL] = 1.0-(-delta-a0)*dt;
        A_CaL[C0_CaL][C1_CaL] = -b0*dt ;
        A_CaL[C0_CaL][C2_CaL] = 0.0;
        A_CaL[C0_CaL][C3_CaL] = 0.0;
        A_CaL[C0_CaL][O_CaL] = 0.0;
        A_CaL[C0_CaL][Ivf_CaL] = 0.0;
        A_CaL[C0_CaL][Ivs_CaL] = 0.0;
        A_CaL[C0_CaL][C0Ca_CaL] = -dt*theta;
        A_CaL[C0_CaL][C1Ca_CaL] = 0.0;
        A_CaL[C0_CaL][C2Ca_CaL] = 0.0;
        A_CaL[C0_CaL][C3Ca_CaL] = 0.0;
        A_CaL[C0_CaL][ICa_CaL] = 0.0;
        A_CaL[C0_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C0_CaL][IvsCa_CaL] = 0.0;

        A_CaL[C1_CaL][C0_CaL] = -a0*dt;
        A_CaL[C1_CaL][C1_CaL] = 1.0-(-delta-b0-a1)*dt  ;
        A_CaL[C1_CaL][C2_CaL] = -b1*dt;
        A_CaL[C1_CaL][C3_CaL] = 0.0;
        A_CaL[C1_CaL][O_CaL] = 0.0;
        A_CaL[C1_CaL][Ivf_CaL] = 0.0;
        A_CaL[C1_CaL][Ivs_CaL] = 0.0;
        A_CaL[C1_CaL][C0Ca_CaL] = 0.0;
        A_CaL[C1_CaL][C1Ca_CaL] = -dt*theta;
        A_CaL[C1_CaL][C2Ca_CaL] = 0.0;
        A_CaL[C1_CaL][C3Ca_CaL] = 0.0;
        A_CaL[C1_CaL][ICa_CaL] = 0.0;
        A_CaL[C1_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C1_CaL][IvsCa_CaL] = 0.0;

        A_CaL[C2_CaL][C0_CaL] = 0.0;
        A_CaL[C2_CaL][C1_CaL] = -a1*dt  ;
        A_CaL[C2_CaL][C2_CaL] =  1.0 -(- delta-b1-a2)*dt;
        A_CaL[C2_CaL][C3_CaL] = -b2*dt;
        A_CaL[C2_CaL][O_CaL] = 0.0;
        A_CaL[C2_CaL][Ivf_CaL] = 0.0;
        A_CaL[C2_CaL][Ivs_CaL] = 0.0;
        A_CaL[C2_CaL][C0Ca_CaL] = 0.0;
        A_CaL[C2_CaL][C1Ca_CaL] = 0.0;
        A_CaL[C2_CaL][C2Ca_CaL] = -dt*theta;
        A_CaL[C2_CaL][C3Ca_CaL]=0.0;
        A_CaL[C2_CaL][ICa_CaL] = 0.0;
        A_CaL[C2_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C2_CaL][IvsCa_CaL] = 0.0;

        A_CaL[C3_CaL][C0_CaL] = 0.0;
        A_CaL[C3_CaL][C1_CaL] = 0.0;
        A_CaL[C3_CaL][C2_CaL] = -a2*dt;
        A_CaL[C3_CaL][C3_CaL] = 1.0-dt*(-ys-yf-delta-b2-a3);
        A_CaL[C3_CaL][O_CaL] = -b3*dt;
        A_CaL[C3_CaL][Ivf_CaL] = -dt*wf;
        A_CaL[C3_CaL][Ivs_CaL] = -dt*ws;
        A_CaL[C3_CaL][C0Ca_CaL] = 0.0;
        A_CaL[C3_CaL][C1Ca_CaL] = 0.0;
        A_CaL[C3_CaL][C2Ca_CaL] = 0.0;
        A_CaL[C3_CaL][C3Ca_CaL]=-dt*theta;
        A_CaL[C3_CaL][ICa_CaL] = 0.0;
        A_CaL[C3_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C3_CaL][IvsCa_CaL] = 0.0;

        A_CaL[O_CaL][C0_CaL] = 0.0;
        A_CaL[O_CaL][C1_CaL] = 0.0;
        A_CaL[O_CaL][C2_CaL] = 0.0;
        A_CaL[O_CaL][C3_CaL] = -a3*dt;
        A_CaL[O_CaL][O_CaL] = 1.0 -dt*(-ps-pf-delta-b3);
        A_CaL[O_CaL][Ivf_CaL] = -dt*lf;
        A_CaL[O_CaL][Ivs_CaL] = -dt*ls;
        A_CaL[O_CaL][C0Ca_CaL] = 0.0;
        A_CaL[O_CaL][C1Ca_CaL] = 0.0;
        A_CaL[O_CaL][C2Ca_CaL] = 0.0;
        A_CaL[O_CaL][C3Ca_CaL] = 0.0;
        A_CaL[O_CaL][ICa_CaL] = -dt*theta;
        A_CaL[O_CaL][IvfCa_CaL] = 0.0;
        A_CaL[O_CaL][IvsCa_CaL] = 0.0;

        A_CaL[Ivf_CaL][C0_CaL] = 0.0;
        A_CaL[Ivf_CaL][C1_CaL] = 0.0;
        A_CaL[Ivf_CaL][C2_CaL] = 0.0;
        A_CaL[Ivf_CaL][C3_CaL] = -dt*yf;
        A_CaL[Ivf_CaL][O_CaL] = -dt*pf;
        A_CaL[Ivf_CaL][Ivf_CaL] = 1.0-dt*(-wfs-wf-lf-delta);
        A_CaL[Ivf_CaL][Ivs_CaL] = -dt*wsf;
        A_CaL[Ivf_CaL][C0Ca_CaL] = 0.0;
        A_CaL[Ivf_CaL][C1Ca_CaL] = 0.0;
        A_CaL[Ivf_CaL][C2Ca_CaL] = 0.0;
        A_CaL[Ivf_CaL][C3Ca_CaL] = 0.0;
        A_CaL[Ivf_CaL][ICa_CaL] = 0.0;
        A_CaL[Ivf_CaL][IvfCa_CaL] = -dt*theta;
        A_CaL[Ivf_CaL][IvsCa_CaL] = 0.0;

        A_CaL[Ivs_CaL][C0_CaL] = 0.0;
        A_CaL[Ivs_CaL][C1_CaL] = 0.0;
        A_CaL[Ivs_CaL][C2_CaL] = 0.0;
        A_CaL[Ivs_CaL][C3_CaL] = -dt*ys;
        A_CaL[Ivs_CaL][O_CaL] = -dt*ps;
        A_CaL[Ivs_CaL][Ivf_CaL] = -dt*wfs;
        A_CaL[Ivs_CaL][Ivs_CaL] = 1.0-dt*(-wsf-ws-ls-delta) ;
        A_CaL[Ivs_CaL][C0Ca_CaL] = 0.0;
        A_CaL[Ivs_CaL][C1Ca_CaL] = 0.0;
        A_CaL[Ivs_CaL][C2Ca_CaL] = 0.0;
        A_CaL[Ivs_CaL][C3Ca_CaL] = 0.0;
        A_CaL[Ivs_CaL][ICa_CaL] = 0.0;
        A_CaL[Ivs_CaL][IvfCa_CaL] = 0.0;
        A_CaL[Ivs_CaL][IvsCa_CaL] = -dt*theta;

        A_CaL[C0Ca_CaL][C0_CaL] = -delta*dt;
        A_CaL[C0Ca_CaL][C1_CaL] = 0.0;
        A_CaL[C0Ca_CaL][C2_CaL] = 0.0;
        A_CaL[C0Ca_CaL][C3_CaL] = 0.0;
        A_CaL[C0Ca_CaL][O_CaL] = 0.0;
        A_CaL[C0Ca_CaL][Ivf_CaL] = 0.0;
        A_CaL[C0Ca_CaL][Ivs_CaL] = 0.0;
        A_CaL[C0Ca_CaL][C0Ca_CaL] = 1.0-dt*(-theta- a0);
        A_CaL[C0Ca_CaL][C1Ca_CaL] = -b0*dt;
        A_CaL[C0Ca_CaL][C2Ca_CaL] = 0.0;
        A_CaL[C0Ca_CaL][C3Ca_CaL] = 0.0;
        A_CaL[C0Ca_CaL][ICa_CaL] = 0.0;
        A_CaL[C0Ca_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C0Ca_CaL][IvsCa_CaL] = 0.0;

        A_CaL[C1Ca_CaL][C0_CaL] = 0.0;
        A_CaL[C1Ca_CaL][C1_CaL] = -delta*dt;
        A_CaL[C1Ca_CaL][C2_CaL] = 0.0;
        A_CaL[C1Ca_CaL][C3_CaL] = 0.0;
        A_CaL[C1Ca_CaL][O_CaL] = 0.0;
        A_CaL[C1Ca_CaL][Ivf_CaL] = 0.0;
        A_CaL[C1Ca_CaL][Ivs_CaL] = 0.0;
        A_CaL[C1Ca_CaL][C0Ca_CaL] = -a0*dt;
        A_CaL[C1Ca_CaL][C1Ca_CaL] = 1.0-dt*(-theta-b0-a1);
        A_CaL[C1Ca_CaL][C2Ca_CaL] = -b1*dt;
        A_CaL[C1Ca_CaL][C3Ca_CaL] = 0.0;
        A_CaL[C1Ca_CaL][ICa_CaL] = 0.0;
        A_CaL[C1Ca_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C1Ca_CaL][IvsCa_CaL] = 0.0;

        A_CaL[C2Ca_CaL][C0_CaL] = 0.0;
        A_CaL[C2Ca_CaL][C1_CaL] = 0.0;
        A_CaL[C2Ca_CaL][C2_CaL] = -delta*dt;
        A_CaL[C2Ca_CaL][C3_CaL] = 0.0;
        A_CaL[C2Ca_CaL][O_CaL] = 0.0;
        A_CaL[C2Ca_CaL][Ivf_CaL] = 0.0;
        A_CaL[C2Ca_CaL][Ivs_CaL] = 0.0;
        A_CaL[C2Ca_CaL][C0Ca_CaL] = 0.0;
        A_CaL[C2Ca_CaL][C1Ca_CaL] = -a1*dt;
        A_CaL[C2Ca_CaL][C2Ca_CaL] = 1.0-dt*(-theta-b1-a2);
        A_CaL[C2Ca_CaL][C3Ca_CaL] = -b2*dt;
        A_CaL[C2Ca_CaL][ICa_CaL] = 0.0;
        A_CaL[C2Ca_CaL][IvfCa_CaL] = 0.0;
        A_CaL[C2Ca_CaL][IvsCa_CaL] = 0.0;

        A_CaL[C3Ca_CaL][C0_CaL] = 0.0;
        A_CaL[C3Ca_CaL][C1_CaL] = 0.0;
        A_CaL[C3Ca_CaL][C2_CaL] = 0.0;
        A_CaL[C3Ca_CaL][C3_CaL] = -delta*dt;
        A_CaL[C3Ca_CaL][O_CaL] = 0.0;
        A_CaL[C3Ca_CaL][Ivf_CaL] = 0.0;
        A_CaL[C3Ca_CaL][Ivs_CaL] = 0.0;
        A_CaL[C3Ca_CaL][C0Ca_CaL] = 0.0;
        A_CaL[C3Ca_CaL][C1Ca_CaL] = 0.0;
        A_CaL[C3Ca_CaL][C2Ca_CaL] = -a2*dt;
        A_CaL[C3Ca_CaL][C3Ca_CaL] = 1.0-dt*(-ys-yf-theta-b2-a3);
        A_CaL[C3Ca_CaL][ICa_CaL] = -b3*dt;
        A_CaL[C3Ca_CaL][IvfCa_CaL] = -dt*wf;
        A_CaL[C3Ca_CaL][IvsCa_CaL] = -dt*ws;

        A_CaL[ICa_CaL][C0_CaL] = 0.0;
        A_CaL[ICa_CaL][C1_CaL] = 0.0;
        A_CaL[ICa_CaL][C2_CaL] = 0.0;
        A_CaL[ICa_CaL][C3_CaL] = 0.0;
        A_CaL[ICa_CaL][O_CaL] = -delta*dt;
        A_CaL[ICa_CaL][Ivf_CaL] = 0.0;
        A_CaL[ICa_CaL][Ivs_CaL] = 0.0;
        A_CaL[ICa_CaL][C0Ca_CaL] = 0.0;
        A_CaL[ICa_CaL][C1Ca_CaL] = 0.0;
        A_CaL[ICa_CaL][C2Ca_CaL] = 0.0;
        A_CaL[ICa_CaL][C3Ca_CaL] = -a3*dt;
        A_CaL[ICa_CaL][ICa_CaL] = 1.0-dt*(-theta-ps-pf-b3);
        A_CaL[ICa_CaL][IvfCa_CaL] = -dt*lf;
        A_CaL[ICa_CaL][IvsCa_CaL] = -dt*ls;

        A_CaL[IvfCa_CaL][C0_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C1_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C2_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C3_CaL] = 0.0;
        A_CaL[IvfCa_CaL][O_CaL] = 0.0;
        A_CaL[IvfCa_CaL][Ivf_CaL] = -delta*dt;
        A_CaL[IvfCa_CaL][Ivs_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C0Ca_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C1Ca_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C2Ca_CaL] = 0.0;
        A_CaL[IvfCa_CaL][C3Ca_CaL] = -dt*yf;
        A_CaL[IvfCa_CaL][ICa_CaL] = -dt*pf;
        A_CaL[IvfCa_CaL][IvfCa_CaL] = 1.0-dt*(-wfs-wf-theta-lf);
        A_CaL[IvfCa_CaL][IvsCa_CaL] = -dt*wsf;

        A_CaL[IvsCa_CaL][C0_CaL] = 0.0;
        A_CaL[IvsCa_CaL][C1_CaL] = 0.0;
        A_CaL[IvsCa_CaL][C2_CaL] = 0.0;
        A_CaL[IvsCa_CaL][C3_CaL] = 0.0;
        A_CaL[IvsCa_CaL][O_CaL] = 0.0;
        A_CaL[IvsCa_CaL][Ivf_CaL] = 0.0;
        A_CaL[IvsCa_CaL][Ivs_CaL] = -delta*dt;
        A_CaL[IvsCa_CaL][C0Ca_CaL] = 0.0;
        A_CaL[IvsCa_CaL][C1Ca_CaL] = 0.0;
        A_CaL[IvsCa_CaL][C2Ca_CaL] = 0.0;
        A_CaL[IvsCa_CaL][C3Ca_CaL] = -dt*ys;
        A_CaL[IvsCa_CaL][ICa_CaL] = -dt*ps;
        A_CaL[IvsCa_CaL][IvfCa_CaL] = -dt*wfs;
        A_CaL[IvsCa_CaL][IvsCa_CaL] = 1.0-dt*(-wsf-ws-theta-ls);

        //Solve using LU decomposition
        ludcmp(A_CaL,n_nr_CaL,indx_CaL,&d_nr_CaL);
        lubksb(A_CaL,n_nr_CaL,indx_CaL,B_CaL);

        //Normalization
        norm_stat_val(n_nr_CaL,B_CaL);

        ICaL=GCaL*B_CaL[O_CaL]*(V-ECa); //pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ---------------------- */
        /* T-type calcium current */
        /* ---------------------- */

        tau_dCaT=1.9508/T_correction_Ca; //msec
        tau_fCaT=(0.38117*(8.6+14.7*exp(-(V+50.0)*(V+50.0)/900.0)))/T_correction_Ca; //msec

        dCaT_inf=1.0/(1.0+exp(-(V+60.5)/5.3));
        fCaT_inf=1.0/(1.0+exp((V+75.5)/4.0));

        ddCaTdt=(dCaT_inf-dCaT)/tau_dCaT;
        dfCaTdt=(fCaT_inf-fCaT)/tau_fCaT;

        dCaT=dCaT+dt*ddCaTdt;
        fCaT=fCaT+dt*dfCaTdt;

        ICaT=GCaT*dCaT*fCaT*(V-ECa); //pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ----------------------------------- */
        /* Voltage dependent potassium current */
        /* ----------------------------------- */

        tau_xKv= 1.0*4.7803/T_correction_K; //msec
        tau_yKv= 1.0*763.7564/T_correction_K; //msec

        xKv_inf=1.0/(1.0+exp(-(V+43.0)/(17.36)));
        yKv_inf=1.0/(1.0+exp((V-44.9)/(12.0096)));

        dxKvdt=(xKv_inf-xKv)/tau_xKv;
        dyKvdt=(yKv_inf-yKv)/tau_yKv;

        xKv=xKv+dt*dxKvdt;
        yKv=yKv+dt*dyKvdt;

        IKv=GKv*xKv*yKv*(V-EK); //pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ------------------------------------------------------- */
        /* Voltage dependent & calcium sensitive potassium current */
        /* ------------------------------------------------------- */

        //Rate equations
        F_BK=6.022*1e23; //mol^(-1), faraday's constant
        R_BK=5.189*1e19; //eVK^(-1) mol^(-1), ideal gas constant
        T_BK=T; //K
        V_m=-1.0*V*1e-3;//convert mV to V, and reverse polarity
        Ca_free=Cai_free*1.0e-3;// convert from mM to M
        time_convert=1.0e-3;//the rate values are in per sec, convert to msec

        qCO=0.73;
        commonval1=1.0*exp((-qCO*F_BK*V_m)/(R_BK*T_BK));//reduce computational cost
        k_C0O0=time_convert*paramBK[0]*commonval1;
        k_C1O1=time_convert*paramBK[1]*commonval1;
        k_C2O2=time_convert*paramBK[2]*commonval1;
        k_C3O3=time_convert*paramBK[3]*commonval1;
        k_C4O4=time_convert*paramBK[4]*commonval1;

        qOC=-0.67;
        commonval2=1.0*exp((-qOC*F_BK*V_m)/(R_BK*T_BK));
        k_O0C0=time_convert*paramBK[5]*commonval2;
        k_O1C1=time_convert*paramBK[6]*commonval2;
        k_O2C2=time_convert*paramBK[7]*commonval2;
        k_O3C3=time_convert*paramBK[8]*commonval2;
        k_O4C4=time_convert*paramBK[9]*commonval2;

        k_ON=1.0e9*paramBK[10]; //M^(-1) s^(-1)
        k_OFF_C=1.0e9*11.0*1.0e-6*paramBK[11];//s^(-1)
        k_OFF_O=1.0e9*1.1*1.0e-6*paramBK[12];//s^(-1)

        k_C0C1=time_convert*4.0*k_ON*Ca_free;
        k_C1C2=time_convert*3.0*k_ON*Ca_free;
        k_C2C3=time_convert*2.0*k_ON*Ca_free;
        k_C3C4=time_convert*k_ON*Ca_free;
        k_C4C3=time_convert*4.0*k_OFF_C;
        k_C3C2=time_convert*3.0*k_OFF_C;
        k_C2C1=time_convert*2.0*k_OFF_C;
        k_C1C0=time_convert*k_OFF_C;

        k_O0O1=time_convert*4.0*k_ON*Ca_free;
        k_O1O2=time_convert*3.0*k_ON*Ca_free;
        k_O2O3=time_convert*2.0*k_ON*Ca_free;
        k_O3O4=time_convert*k_ON*Ca_free;
        k_O4O3=time_convert*4.0*k_OFF_O;
        k_O3O2=time_convert*3.0*k_OFF_O;
        k_O2O1=time_convert*2.0*k_OFF_O;
        k_O1O0=time_convert*k_OFF_O;

        //Equation system rearranged into Ax=b format, following backward Euler implementation
        //Assign the row vectors of the A matrix
        A_BK[O0_BK][O0_BK] = 1.0 -dt*(-k_O0O1-k_O0C0);
        A_BK[O0_BK][O1_BK] = -dt*k_O1O0;
        A_BK[O0_BK][O2_BK] = 0.0;
        A_BK[O0_BK][O3_BK] = 0.0;
        A_BK[O0_BK][O4_BK] = 0.0;
        A_BK[O0_BK][C0_BK] = -dt*k_C0O0;
        A_BK[O0_BK][C1_BK] = 0.0;
        A_BK[O0_BK][C2_BK] = 0.0;
        A_BK[O0_BK][C3_BK] = 0.0;
        A_BK[O0_BK][C4_BK] = 0.0;

        A_BK[O1_BK][O0_BK] = -dt*k_O0O1  ;
        A_BK[O1_BK][O1_BK] = 1.0 - dt*(-k_O1O2-k_O1O0-k_O1C1);
        A_BK[O1_BK][O2_BK] = -dt*k_O2O1 ;
        A_BK[O1_BK][O3_BK] = 0.0;
        A_BK[O1_BK][O4_BK] = 0.0;
        A_BK[O1_BK][C0_BK] = 0.0;
        A_BK[O1_BK][C1_BK] = -dt*k_C1O1 ;
        A_BK[O1_BK][C2_BK] = 0.0;
        A_BK[O1_BK][C3_BK] = 0.0;
        A_BK[O1_BK][C4_BK] = 0.0;

        A_BK[O2_BK][O0_BK] = 0.0;
        A_BK[O2_BK][O1_BK] = -dt*k_O1O2  ;
        A_BK[O2_BK][O2_BK] = 1.0-dt*(-k_O2O3-k_O2O1-k_O2C2);
        A_BK[O2_BK][O3_BK] = -dt*k_O3O2;
        A_BK[O2_BK][O4_BK] = 0.0;
        A_BK[O2_BK][C0_BK] = 0.0;
        A_BK[O2_BK][C1_BK] = 0.0;
        A_BK[O2_BK][C2_BK] = -dt*k_C2O2 ;
        A_BK[O2_BK][C3_BK] = 0.0;
        A_BK[O2_BK][C4_BK] = 0.0;

        A_BK[O3_BK][O0_BK] = 0.0;
        A_BK[O3_BK][O1_BK] = 0.0;
        A_BK[O3_BK][O2_BK] = -dt*k_O2O3  ;
        A_BK[O3_BK][O3_BK] = 1.0-dt*(-k_O3O4-k_O3O2-k_O3C3);
        A_BK[O3_BK][O4_BK] = -dt*k_O4O3   ;
        A_BK[O3_BK][C0_BK] = 0.0;
        A_BK[O3_BK][C1_BK] = 0.0;
        A_BK[O3_BK][C2_BK] = 0.0;
        A_BK[O3_BK][C3_BK] = -dt*k_C3O3 ;
        A_BK[O3_BK][C4_BK] = 0.0;

        A_BK[O4_BK][O0_BK] = 0.0;
        A_BK[O4_BK][O1_BK] = 0.0;
        A_BK[O4_BK][O2_BK] = 0.0;
        A_BK[O4_BK][O3_BK] = -dt*k_O3O4   ;
        A_BK[O4_BK][O4_BK] = 1.0-dt*(-k_O4O3-k_O4C4);
        A_BK[O4_BK][C0_BK] = 0.0;
        A_BK[O4_BK][C1_BK] = 0.0;
        A_BK[O4_BK][C2_BK] = 0.0;
        A_BK[O4_BK][C3_BK] = 0.0;
        A_BK[O4_BK][C4_BK] = -dt*k_C4O4;

        A_BK[C0_BK][O0_BK] = -dt*k_O0C0;
        A_BK[C0_BK][O1_BK] = 0.0;
        A_BK[C0_BK][O2_BK] = 0.0;
        A_BK[C0_BK][O3_BK] = 0.0;
        A_BK[C0_BK][O4_BK] = 0.0;
        A_BK[C0_BK][C0_BK] = 1.0-dt*(-k_C0O0-k_C0C1);
        A_BK[C0_BK][C1_BK] = -dt*k_C1C0;
        A_BK[C0_BK][C2_BK] = 0.0;
        A_BK[C0_BK][C3_BK] = 0.0;
        A_BK[C0_BK][C4_BK] = 0.0;

        A_BK[C1_BK][O0_BK] = 0.0;
        A_BK[C1_BK][O1_BK] = -dt*k_O1C1;
        A_BK[C1_BK][O2_BK] = 0.0;
        A_BK[C1_BK][O3_BK] = 0.0;
        A_BK[C1_BK][O4_BK] = 0.0;
        A_BK[C1_BK][C0_BK] = -dt*k_C0C1;
        A_BK[C1_BK][C1_BK] = 1.0-dt*(-k_C1O1-k_C1C2-k_C1C0);
        A_BK[C1_BK][C2_BK] = -dt*k_C2C1 ;
        A_BK[C1_BK][C3_BK] = 0.0;
        A_BK[C1_BK][C4_BK] = 0.0;

        A_BK[C2_BK][O0_BK] = 0.0;
        A_BK[C2_BK][O1_BK] = 0.0;
        A_BK[C2_BK][O2_BK] = -dt*k_O2C2;
        A_BK[C2_BK][O3_BK] = 0.0;
        A_BK[C2_BK][O4_BK] = 0.0;
        A_BK[C2_BK][C0_BK] = 0.0;
        A_BK[C2_BK][C1_BK] = -dt*k_C1C2;
        A_BK[C2_BK][C2_BK] = 1.0-dt*(-k_C2O2-k_C2C3-k_C2C1);
        A_BK[C2_BK][C3_BK] = -dt*k_C3C2;
        A_BK[C2_BK][C4_BK] = 0.0;

        A_BK[C3_BK][O0_BK] = 0.0;
        A_BK[C3_BK][O1_BK] = 0.0;
        A_BK[C3_BK][O2_BK] = 0.0;
        A_BK[C3_BK][O3_BK] =-dt*k_O3C3;
        A_BK[C3_BK][O4_BK] = 0.0;
        A_BK[C3_BK][C0_BK] = 0.0;
        A_BK[C3_BK][C1_BK] = 0.0;
        A_BK[C3_BK][C2_BK] = -dt*k_C2C3 ;
        A_BK[C3_BK][C3_BK] = 1.0-dt*(-k_C3O3-k_C3C4-k_C3C2);
        A_BK[C3_BK][C4_BK] = -dt*k_C4C3  ;

        A_BK[C4_BK][O0_BK] = 0.0;
        A_BK[C4_BK][O1_BK] = 0.0;
        A_BK[C4_BK][O2_BK] = 0.0;
        A_BK[C4_BK][O3_BK] = 0.0;
        A_BK[C4_BK][O4_BK] = -dt*k_O4C4;
        A_BK[C4_BK][C0_BK] = 0.0;
        A_BK[C4_BK][C1_BK] = 0.0;
        A_BK[C4_BK][C2_BK] = 0.0;
        A_BK[C4_BK][C3_BK] = -dt*k_C3C4;
        A_BK[C4_BK][C4_BK] = 1.0-dt*(-k_C4O4-k_C4C3);

        //Solve using LU decomposition
        ludcmp(A_BK,n_nr_BK,indx_BK,&d_nr_BK);
        lubksb(A_BK,n_nr_BK,indx_BK,B_BK);

        //Normalization
    norm_stat_val(n_nr_BK,B_BK);

        IBK=GBK*B_BK[O4_BK]*(V-EK); //pA. O4 is defined to the open state that conducts K ions.

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */

        /* -------------- */
        /* Sodium current */
        /* -------------- */

        //Rate equations for the markov model. all exp
        k_O_I1=T_correction_Na*parval[22]*exp((parval [0]+parval [1]*V));//inact
        k_I1_I2=T_correction_Na*parval[23]*exp((parval [2]+parval [3]*V));//inact
        k_C3_C2=T_correction_Na*parval[24]*exp((parval [4]+parval [5]*V));
        k_C2_C1=T_correction_Na*parval[25]*exp((parval [6]+parval [7]*V));

        k_C1_O=T_correction_Na*parval[26]*exp((parval [8]+parval [9]*V));
        k_I2_I1=T_correction_Na*parval[27]*exp((parval [10]+parval [11]*V));//inact
        k_C2_C3=T_correction_Na*parval[28]*exp((parval [12]+parval [13]*V));
        k_C1_C2=T_correction_Na*parval[29]*exp((parval [14]+parval [15]*V));

        k_O_C1=T_correction_Na*parval[30]*exp((parval [16]+parval [17]*V));
        k_I1_C1=T_correction_Na*parval[31]*exp((parval [18]+parval [19]*V));//inact
        k_C1_I1=T_correction_Na*parval[32]*exp((parval [20]+parval [21]*V));//inact
        k_I1_O=T_correction_Na*parval[33]*exp((parval [34]+parval [35]*V));//inact

        //Equation system rearranged into Ax=b format, following backward Euler implementation
        //Assign the row vectors of the A matrix

        A[ONa][ONa] = 1.0- dt*(- (k_O_I1 + k_O_C1));
        A[ONa][C1] = -k_C1_O*dt;
        A[ONa][C2] = 0.0;
        A[ONa][C3] = 0.0;
        A[ONa][I1] = -k_I1_O*dt;
        A[ONa][I2] = 0.0;

        A[C1][ONa] = -k_O_C1*dt;
        A[C1][C1] = 1.0 - dt*( - (k_C1_O + k_C1_I1 + k_C1_C2));
        A[C1][C2] = -dt*k_C2_C1;
        A[C1][C3] = 0.0;
        A[C1][I1] = -dt*k_I1_C1;
        A[C1][I2] = 0.0;

        A[C2][ONa] = 0;
        A[C2][C1] = -dt*k_C1_C2;
        A[C2][C2] = 1.0-dt*(- (k_C2_C1 + k_C2_C3));
        A[C2][C3] = -dt*k_C3_C2;
        A[C2][I1] = 0.0;
        A[C2][I2] = 0.0;

        A[C3][ONa] = 0.0;
        A[C3][C1] = 0.0;
        A[C3][C2] = -dt*k_C2_C3;
        A[C3][C3] = 1.0+dt*k_C3_C2;
        A[C3][I1] = 0.0;
        A[C3][I2] = 0.0;

        A[I1][ONa] = -dt*k_O_I1;
        A[I1][C1] = -dt*k_C1_I1;
        A[I1][C2] = 0.0;
        A[I1][C3] = 0.0;
        A[I1][I1] = 1.0-dt*( - (k_I1_C1 + k_I1_I2 + k_I1_O));
        A[I1][I2] = -dt*k_I2_I1;

        A[I2][ONa] = 0.0;
        A[I2][C1] = 0.0;
        A[I2][C2] = 0.0;
        A[I2][C3] = 0.0;
        A[I2][I1] = -dt*k_I1_I2;
        A[I2][I2] = 1.0+dt*k_I2_I1;

        //Solve using LU decomposition
        ludcmp(A,n_nr,indx,&d_nr);
        lubksb(A,n_nr,indx,B);

    //Normalization
    norm_stat_val(n_nr,B);

        INa=GNa*B[ONa]*(V-ENa); //pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* --------------------- */
        /* Sodium potassium pump */
        /* --------------------- */

        KmK=1.0;//mM
        KmNa1=40.0;//mM
        numINaK=Ko*Nai;
        denomINaK=(Ko+KmK)*(Nai+KmNa1)*(1+0.1245*exp(-0.1*V*F/(R*T))+0.0353*exp(-V*F/(R*T)));
        INaK=PNaK*numINaK/denomINaK;//pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ------------------------ */
        /* Sodium calcium exchanger */
        /* ------------------------ */

    gamma=0.35;
        KmNa2=87.5;//mM
        KmCa=1.38;//mM
    Ksat=0.1;
        numINCX=exp(0.35*V*F/(R*T))*Nai*Nai*Nai*Cao-2.5*exp(-0.65*V*F/(R*T))*Nao*Nao*Nao*Cai_free;
        denomINCX=(KmNa2*KmNa2*KmNa2+Nao*Nao*Nao)*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*V*F/(R*T)));
        INCX=KNCX*(numINCX/denomINCX); //pA

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* ---------------------- */
        /* Cytosolic Ca buffering */
        /* ---------------------- */
    //Assumption is that buffering reactions reached equilibrium
    //Dissociation constants based derivations
    theta_1_cabuff = (hill_coeff_CRT*CRT_max*(pow(Cai_free, (hill_coeff_CRT - 1)))*
              (K_CRT))/(pow((pow(Cai_free,hill_coeff_CRT) + K_CRT),2));

    theta_2_cabuff = (hill_coeff_CaM*CaM_max*(pow(Cai_free, (hill_coeff_CaM - 1)))*
              (K_CaM))/(pow((pow(Cai_free,hill_coeff_CaM) + K_CaM),2));

    theta_cabuff = 1.0 + theta_1_cabuff + theta_2_cabuff;

    /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* -------------- */
        /* Concentrations */
        /* -------------- */

    //1M=1mol/liter; 1mm^3=1e-6liter;1e-3m^3=1liter
    //ionic currents are in pA

    //Track Cai
    dCaidt=(1.0/(z*F*Vcell))*(-(ICaL+ICaT-2.0*INCX));//Total Cai
        Cai=Cai+(dt*dCaidt)*1.0e-9;//mM, total Cai
        Cai_free=Cai_free+(dt*dCaidt)*1.0e-9/theta_cabuff;//mM, free Cai
    ECa=(R*T/(2.0*F))*log(Cao/Cai_free); //mV

        //Track Nai
    dNaidt=(1.0/(F*Vcell))*-(INa+3.0*INaK+3.0*INCX+gNSNa*(V-ENa));
        Nai=Nai+(dt*dNaidt)*1.0e-9;//mM
    ENa=(R*T/F)*log(Nao/Nai); //mV

    //Track Ki
    dKidt=(1.0/(F*Vcell))*-(IKv+IBK-Istim-2.0*INaK+gNSK*(V-EK));
    Ki=Ki+(dt*dKidt)*1.0e-9;//mM
    EK=(R*T/F)*log(Ko/Ki); //mV

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* --------------------- */
        /* JSMC membrane voltage */
        /* --------------------- */

        INSCC=0.0;//Switch off INSCC
        dVdt=-(1.0/Cm)*(ICaL+ICaT+IKv+IBK+INa+INCX+INaK-Istim+INS_K+INS_Na);
        V=V+dt*dVdt; //mV

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* --------------------------- */
        /* Print the results into file */
        /* --------------------------- */

        printfidx=(i%printfrequency);
        if ((printfidx==0)&&(i>=((100*Slow_wave_period)+Init_equilibration)/dt))
        fprintf (pResults, "%f\t%f\t\%f\t%f\t%f\t%f\t%f\t%f\t\%f\t%f\t%f\t%f\t%f\t%f\t%f\n",1e-3*(t-0),
              V,ICaL,ICaT,IKv,IBK,INa,INCX,INaK,Istim,VICC,Cai,Cai_free,ENa,EK);

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */


        /* --------------------------- */
        /* Print the results on screen */
        /* --------------------------- */

        printsidx=(i%(20*(printfrequency/10)));
        if ((printsidx==0)&&(i>=(Init_equilibration)/dt))
        printf( "Iter %d. Max Iter %f. V %f. \n", i, ((Num_slow_waves*Slow_wave_period)+Init_equilibration)/dt, V);

        /* --------------------------------------------------------------- */
        /* --------------------------------------------------------------- */

        i++;

}
/* ---------end of time loop------------ */


//Free the NR variables of the Markov models
//For CaL model
free_ivector(indx_CaL,1,n_nr_CaL);
free_dvector(B_CaL,1,n_nr_CaL);
free_dmatrix(A_CaL,1,n_nr_CaL,1,n_nr_CaL);
//For BK model
free_ivector(indx_BK,1,n_nr_BK);
free_dvector(B_BK,1,n_nr_BK);
free_dmatrix(A_BK,1,n_nr_BK,1,n_nr_BK);
//For Na model
free_ivector(indx,1,n_nr);
free_dvector(B,1,n_nr);
free_dmatrix(A,1,n_nr,1,n_nr);

fclose(pResults); //close the results file

return 0;
}

void norm_stat_val(int n, double b[])
{
        int i;
        double sum=0.0;
    for (i=1;i<=n;i++) sum=sum+b[i];
        for (i=1;i<=n;i++) b[i]=b[i]/sum;
}
