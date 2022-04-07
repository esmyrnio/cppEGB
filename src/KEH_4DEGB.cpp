#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <string.h>

#include "../include/savgol.hpp"


#define SMAX 1                /* maximum value of s-coordinate */  
#define DS (SMAX/(SDIV-1.0))  

#define KAPPA 1.346790806509621e+13  /* square of length scale = 1e-15*C*C/G */
#define KSCALE 1.112668301525780e-36 /* KAPPA*G/(C*C*C*C) */  
#define RMIN 1.0e-15

#define SDIV 1301
#define RDIV 901 

#define G 6.6732e-8
#define C 2.9979e+10
#define MSUN 1.989e+33
#define PI 3.1415926535

char eos_name[80];

int MAX_ITER, iter;

double Length = G*MSUN/pow(C,2),
Time = Length/C,
Density = MSUN/pow(Length,3),
central_density,
central_pressure,
central_enthalpy,
log_e_tab[201],               /* rho points in tabulated EOS */
log_p_tab[201],               /* p points in tabulated EOS */
log_h_tab[201],              /* h points in EOS file */
log_n0_tab[201],         /* number density in EOS file */  
log_rho_tab[201],
e_tab[201],
p_tab[201],
h_tab[201],
n0_tab[201],
rho_tab[201],
scalar_charge,
r_final=0.0,               
m_final=0.0,
r_is_final=0.0,
r_gp[RDIV+1],
r_is_gp[RDIV+1],
m_gp[RDIV+1],
e_d_gp[RDIV+1],
k_rescale,
re2,
rel_error,
s_e=0.5,
enthalpy_min = 1.0/(C*C),
Mass,
Radius,
RLX_FAC,
TOL;  


double coupling;

double e_surface=7.8*C*C*KSCALE,
        p_surface=1.01e8*KSCALE;

double s_gp[SDIV],
	mu[SDIV],
	nu[SDIV],
	phi[SDIV],
	pressure[SDIV],
	energy_den[SDIV],
	enthalpy[SDIV],
	source_mu[SDIV],
	source_nu[SDIV],
	source_phi[SDIV],
	mu_hat[SDIV],
	nu_hat[SDIV],
	phi_hat[SDIV],
	mu_guess[SDIV],
	nu_guess[SDIV],
	phi_guess[SDIV],
	pressure_guess[SDIV],
	energy_den_guess[SDIV],
	enthalpy_guess[SDIV],
  mu_eq,
  nu_eq,
  nu_hat_eq;

std::vector<double> dmu(SDIV, 0.0);
std::vector<double> dnu(SDIV, 0.0);
std::vector<double> dphi(SDIV, 0.0);
std::vector<double> d2mu(SDIV, 0.0);
std::vector<double> d2nu(SDIV, 0.0);
std::vector<double> d2phi(SDIV, 0.0);
std::vector<double> integral_mu1(SDIV, 0.0);
std::vector<double> integral_mu2(SDIV, 0.0);
std::vector<double> integral_nu1(SDIV, 0.0);
std::vector<double> integral_nu2(SDIV, 0.0);
std::vector<double> integral_phi(SDIV, 0.0);
                  
std::vector<double> dmu_filtered,
                    dnu_filtered,
                    dphi_filtered,
                    d2mu_filtered,
                    d2nu_filtered,
                    d2phi_filtered;
  
double r_e,
	r_e_guess;

double nu_gp[RDIV+1],
	lambda_gp[RDIV+1];

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
void make_grid(void)                                    
{ 
  int s;          /* counter */
   
 for(s=1;s<=SDIV;s++) {
    s_gp[s-1]=(1.0*s-1.0)/(SDIV-1);

  }

}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
double firstOrderDerivative(int i, double f[], double crit){

  if (i==0)
  {
    return crit*(f[i+1]-f[i])/DS;
  }
  else if (i==SDIV-1) 
  {
    return crit*(f[i]-f[i-1])/DS;
  }
  else
  {
    return crit*(f[i+1]-f[i-1])/(2*DS);
  }
}


/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
double trapz(std::vector<double> ar, double size, double step){

  double res;
  double sum = 0;

  for(double x:ar)
    sum+=x;

  res = step*(sum+(ar[0]+ar[size-1])/2);

  return res;
}
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double source_nu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				double e_aw, double p_aw, double r_e, double a){

	double c = 4*PI;
    
  double p = p_aw;
  double eps = e_aw;
  
  double emu = exp(2*mu_t);
  double memu = exp(-2*mu_t);
  
  double r = sqrt(r_e)*s/(1-s);
  double r2 = r_e*pow(s/(1-s),2);
  
  double ss = pow(1-s,2)/sqrt(r_e);
  double dn = ss*dnds;
  double dm = ss*dmds;
  double dphi = ss*dfds;
  
  double ddphi = pow(1-s,4)*dfds2/r_e - dfds*2*pow(1-s,3)/r_e;
  double ddm =   pow(1-s,4)*dmds2/r_e - dmds*2*pow(1-s,3)/r_e;
  
  double term = 2*dn*r;

  if(i!=0){

  	return r*(emu*(c*emu*r2*(3*p+eps)-r*dn*(2+r*(dm+dn))+a*memu*
              (dphi*(2*r2*pow(dm,3)-2*r2*pow(dm,2)*(dn+dphi)-2*r*pow(dn,2)
                *(-2+r*dphi)-r*(pow(dphi,2)*(-4+r*dphi)+4*ddm)+2*dm
              *(-2+r*(2*r*pow(dn,2)+dn*(4-r*dphi)+dphi*(-2+r*dphi)-2*r*ddm))
                +2*dn*(2+r*(dphi*(-4+r*dphi)+2*r*ddm)))+2*r*(-r*pow(dm,2)+2*dm*
                (-1+r*dn)+r*pow(dphi,2)+dn*(2-2*r*dphi))*ddphi)))
              /(emu*r+2*a*dphi*(-2-2*r*dm+r*dphi)) + term ;
      }
  else{
  	return 0;
  	}
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double source_mu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				double e_aw, double p_aw, double r_e, double a){

	double c = 8*PI;
    
  double p = p_aw;
  double eps = e_aw;
  
  double emu = exp(2*mu_t);

  double r = sqrt(r_e)*s/(1-s);
  double r2 = r_e*pow(s/(1-s),2);
  
  double ss = pow(1-s,2)/sqrt(r_e);
  double dn = ss*dnds;
  double dm = ss*dmds;
  double dphi = ss*dfds;
  
  double ddphi = pow(1-s,4)*dfds2/r_e - dfds*2*pow(1-s,3)/r_e;
  
  double term = 2*dm*r;

  if(i!=0){

      return -r*((c*pow(emu,2)*r2*eps+emu*r*dm*(4+r*dm)+a*(2*dm-dphi)
       *dphi*(-4+r2*(2*pow(dm,2)-2*dm*dphi+pow(dphi,2)))
     -4*a*r*(dm-dphi)*(2+r*dm-r*dphi)*ddphi)/
       (2*(emu*r+2*a*dphi*(-2-2*r*dm+r*dphi)))) + term ;
  }
  else{
  	return 0;
  }
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double source_phi_function(int i, double s, double dmds, double r_e){

	double ss = pow(1-s,2)/sqrt(r_e);
	double dm = ss*dmds;

	return dm;
}
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
int n_tab,                           /* number of tabulated EOS points */
    n_nearest=1,                     /* nearest grid point, used in interp. */ 
    print_option;                    /* select print out */ 

void load_eos( char eos_file[])
{
 int i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */    
      
 FILE *f_eos;              /* pointer to eos_file */
 

    /* OPEN FILE TO READ */

    char eos[100];
    strcat(eos,"EoS/");  /* add path */
    strcat(eos,eos_file);
	
    if((f_eos=fopen(eos,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",eos_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",&n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS (EXCEPT N0) */
 
    for(i=1;i<=n_tab;i++) {  
      /*  input line for _four_column equation of state file */
       fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ; 

/*  input line for _five_ column equation of state file         
    fscanf(f_eos,"%lf %lf %lf %lf %*lf\n",&rho,&p,&h,&n0) ; */

       log_e_tab[i]=log10(rho*C*C*KSCALE);       /* multiply by C^2 to get */ 
       log_p_tab[i]=log10(p*KSCALE);             /* energy density. */
       log_h_tab[i]=log10(h/(C*C));        
       log_n0_tab[i]=log10(n0);                  /* STILL IN CGS ! */

    }
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

void hunt(double xx[], int n, double x, int *jlo)
{ 
  int jm,jhi,inc,ascnd;

  ascnd=(xx[SDIV] > xx[1]);
  if (*jlo <= 0 || *jlo > n) {
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    if (x >= xx[*jlo] == ascnd) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while (x >= xx[jhi] == ascnd) {
        *jlo=jhi;
        inc += inc;
        jhi=(*jlo)+inc;
        if (jhi > n) {
          jhi=n+1;
          break;
        }
      }
    } else {
      if (*jlo == 1) {
        *jlo=0;
        return;
      }
      jhi=(*jlo);
      *jlo -= 1;
      while (x < xx[*jlo] == ascnd) {
        jhi=(*jlo);
        inc += inc;
        *jlo=jhi-inc;
        if (*jlo < 1) {
          *jlo=0;
          break;
        }
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if (x > xx[jm] == ascnd)
      *jlo=jm;
    else
      jhi=jm;
  }
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double interp(double xp[], double yp[], int np ,double xb)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,&n_nearest);

 k=std::min(std::max(n_nearest-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) xb += 1e-12;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double e_at_p(double pp)
{
  if(pp<=0) return 0;
  else
    return pow(10.0,interp(log_p_tab,log_e_tab,n_tab,log10(pp)));
}
/*******************************************************************/
double e_at_rho(double rhorho)
{
   return pow(10.0,interp(log_rho_tab,log_e_tab,n_tab,log10(rhorho))); 
}
/*******************************************************************/
double p_at_e(double ee)
{
  return pow(10.0,interp(log_e_tab,log_p_tab,n_tab,log10(ee)));
}
/*******************************************************************/
double p_at_rho(double rhorho)
{
  return pow(10.0,interp(log_rho_tab,log_p_tab,n_tab,log10(rhorho)));
}
/*******************************************************************/
double p_at_h(double hh)
{ 
  if(hh<0) return 0;
    
  else
      return pow(10.0,interp(log_h_tab,log_p_tab,n_tab,log10(hh)));
}
/*******************************************************************/
double h_at_p(double pp)
{ 
  return pow(10.0,interp(log_p_tab,log_h_tab,n_tab,log10(pp)));
}
/*******************************************************************/
double n0_at_e(double ee)
{ 
  return pow(10.0,interp(log_e_tab,log_n0_tab,n_tab,log10(ee)));
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

void make_center(double central_e){
  n_nearest=n_tab/2;
  central_pressure = p_at_e(central_e);
  central_enthalpy = h_at_p(central_pressure);
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
template <typename T>
std::vector<T> slice(std::vector<T> const& v,
                  int X, int Y)
{
    auto first = v.begin() + X;
    auto last = v.begin() + Y + 1;
 
    std::vector<T> vector(first, last);
 
    return vector;
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double dm_dr_is(double r_is, double r, double m, double p)
{
 double dmdr,
        e_d;

 if(p<p_surface) 
    e_d=0.0;
 else 
    e_d=e_at_p(p);
 
 if(r_is<RMIN) 
    dmdr=4.0*PI*central_density*r*r*(1.0+4.0*PI*central_density*r*r/3.0);
 else
    dmdr=4.0*PI*e_d*r*r*r*sqrt(1.0-2.0*m/r)/r_is;
 
return dmdr;
}
 
/**************************************************************************/
double dp_dr_is(double r_is, double r, double m, double p)

{ double dpdr,
         e_d; 

  if(p<p_surface) e_d=0.0;
  else        
   e_d=e_at_p(p);
  
  if(r_is<RMIN) dpdr = -4.0*PI*(central_density+p)*(central_density+3.0*p)*r*(1.0
                     +4.0*central_density*r*r/3.0)/3.0;

  else 
   dpdr = -(e_d+p)*(m+4.0*PI*r*r*r*p)/(r*r_is*sqrt(1.0-2.0*m/r));

 return dpdr;
}

/**************************************************************************/
double dr_dr_is(double r_is, double r, double m)
{
 double drdris;

 if(r_is<RMIN) drdris=1.0;
  else
   drdris=(r/r_is)*sqrt(1.0-2.0*m/r);

 return drdris;
}
/**************************************************************************/
/**************************************************************************/

void integrate(int i_check)
{
  int i=2;

  double r,                           /* radius */
         r_is,                        /* isotropic radial coordinate */
         m,                           /* mass   */
         p,                           /* pressure */
         e_d,                         /* density */
         r_is_est,                    /* estimate on final isotr. radius */ 
         dr_is_save,                  /* r_is saving interval */  
         r_is_check,                  /*                      */    
         nu_s,
         hh,
         h,                           /* stepsize during integration */
         a1,a2,a3,a4,b1,b2,b3,b4,     /* coeff. in Runge-Kutta equations */
         c1,c2,c3,c4,
         rho0;   

    if(i_check==1) {
        r_is_est=1.5e6/sqrt(KAPPA);

        h=r_is_est/100;     
    }
    else {
          r_is_est=r_is_final;
          h=r_is_est/10000;     
   }

    dr_is_save=r_is_final/RDIV;
    r_is_check=dr_is_save;
    r_is=0.0;                            /* initial isotropic radius */
    r=0.0;                               /* initial radius */
    m=0.0;                               /* initial mass */ 
    p=central_pressure;                  /* initial pressure */ 

    r_is_gp[1]=0.0;
    r_gp[1]=0.0;
    m_gp[1]=0.0;
    lambda_gp[1]=0.0;
    e_d_gp[1]=central_density; 

    while (p>=p_surface) { 

      e_d=e_at_p(p);

     if((i_check==3) && (r_is>r_is_check) && (i<=RDIV)) {
      r_is_gp[i]=r_is;
      r_gp[i]=r;
      m_gp[i]=m;
      e_d_gp[i]=e_d; 
      i++;   
      r_is_check += dr_is_save;
     }    
       
     r_is_final=r_is;
     r_final=r;
     m_final=m;
 
     a1=dr_dr_is(r_is,r,m);
     b1=dm_dr_is(r_is,r,m,p);
     c1=dp_dr_is(r_is,r,m,p);

 
     a2=dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0);
     b2=dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0);
     c2=dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0);

     a3=dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0);
     b3=dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0);
     c3=dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0);

     a4=dr_dr_is(r_is+h, r+h*a3, m+h*b3);
     b4=dm_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3);
     c4=dp_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3);

     
     r += (h/6.0)*(a1+2*a2+2*a3+a4);
     m += (h/6.0)*(b1+2*b2+2*b3+b4);
     p += (h/6.0)*(c1+2*c2+2*c3+c4);

     r_is += h;
    }

    r_is_gp[RDIV]=r_is_final;
    r_gp[RDIV]=r_final;
    m_gp[RDIV]=m_final;

/* Rescale r_is and compute lambda */

   if(i_check==3) {
      k_rescale=0.5*(r_final/r_is_final)*(1.0-m_final/r_final+
                sqrt(1.0-2.0*m_final/r_final));
 

      r_is_final *= k_rescale;

      nu_s = log((1.0-m_final/(2.0*r_is_final))/(1.0+m_final/
               (2.0*r_is_final)));

      for(i=1;i<=RDIV;i++) {
         r_is_gp[i] *= k_rescale;
 
         if(i==1) lambda_gp[1]= log(1.0/k_rescale);
           else
               lambda_gp[i]=log(r_gp[i]/r_is_gp[i]); 

         if(e_d_gp[i]<e_surface) 
 
           hh=0.0;
 
         else { 
                  p=p_at_e(e_d_gp[i]);
                  hh=h_at_p(p);
              }
         nu_gp[i]=nu_s-hh;
      }
      nu_gp[RDIV]=nu_s;

   }
}

void guess(void)
{
 int i,
     s,
     m,
     s_temp;

 double r_is_s,
        lambda_s,
        nu_s,
        gama_eq,
        rho_eq;

 n_nearest=n_tab/2;

 /* SOLVE THE OPPENHEIMER-VOLKOV EQUATIONS USING A RUNGE-KUTTA METHOD
    The functions integrate solves the equations using the r_is coordinate */
 integrate(1);
 integrate(2);
 integrate(3);

 if(SMAX==1.0) s_temp=SDIV-1;
 else 
     s_temp=SDIV;


 n_nearest=RDIV/2;

 for(s=0;s<=s_temp;s++) {
    r_is_s=r_is_final*(s_gp[s]/(1.0-s_gp[s]));
    /* CONVERT THE SPHERICAL SOLUTION TO THE "s" "radial" coordinate */
    if(r_is_s<r_is_final) {
      lambda_s=interp(r_is_gp,lambda_gp,RDIV,r_is_s);
      nu_s=interp(r_is_gp,nu_gp,RDIV,r_is_s);
    }
    else {
      lambda_s=2.0*log(1.0+m_final/(2.0*r_is_s));
      nu_s=log((1.0-m_final/(2.0*r_is_s))/(1.0+m_final/(2*r_is_s)));
    }
    nu_guess[s]=nu_s;
    mu_guess[s]=lambda_s;
    phi_guess[s]=lambda_s;
 }

   mu_eq=mu_guess[int(SDIV/2)];

   r_e_guess= r_final*exp(-mu_eq);
}


void main_iteration(void){

  iter=0;

  std::size_t i;

  while(1>0){

    for (i = 0; i < SDIV; ++i)
    {
      mu_hat[i] = mu_guess[i]/pow(r_e_guess,2);
      nu_hat[i] = nu_guess[i]/pow(r_e_guess,2);
      phi_hat[i] = phi_guess[i]/pow(r_e_guess,2);
    }

    nu_hat_eq = nu_hat[int(SDIV/2)];

    r_e = central_enthalpy/(nu_hat_eq-nu_hat[0]);

    n_nearest=n_tab/2;

    for (i = 0; i < SDIV; ++i)
    {
      enthalpy[i]=enthalpy_min + r_e*(nu_hat_eq-nu_hat[i]);
      pressure[i] = p_at_h(enthalpy[i]);
      energy_den[i] = e_at_p(pressure[i]);

      dmu[i] = firstOrderDerivative(i,mu_hat,r_e);
      dnu[i] = firstOrderDerivative(i,nu_hat,r_e);
      dphi[i] = firstOrderDerivative(i,phi_hat,r_e);
    }

    dmu_filtered = sg_smooth(dmu,31,3);
    dnu_filtered = sg_smooth(dnu,31,3);
    dphi_filtered = sg_smooth(dphi,31,3);

    double* dmu_filtered_ar = &dmu_filtered[0];
    double* dnu_filtered_ar = &dnu_filtered[0];
    double* dphi_filtered_ar = &dphi_filtered[0];

    for (i = 0; i < SDIV; ++i)
    {

      d2mu[i] = firstOrderDerivative(i,dmu_filtered_ar,1);
      d2nu[i] = firstOrderDerivative(i,dnu_filtered_ar,1);
      d2phi[i] = firstOrderDerivative(i,dphi_filtered_ar,1);

    }

    d2mu_filtered = sg_smooth(d2mu,31,3);
    d2nu_filtered = sg_smooth(d2nu,31,3);
    d2phi_filtered = sg_smooth(d2phi,31,3);

    double* d2mu_filtered_ar = &d2mu_filtered[0];
    double* d2nu_filtered_ar = &d2nu_filtered[0];
    double* d2phi_filtered_ar = &d2phi_filtered[0];


    for (i = 0; i < SDIV; ++i)
    {
      source_nu[i] = source_nu_function(i, s_gp[i], r_e*mu_hat[i], r_e*nu_hat[i], r_e*phi_hat[i], dnu_filtered_ar[i], d2nu_filtered_ar[i],
                     dmu_filtered_ar[i], d2mu_filtered_ar[i], dphi_filtered_ar[i], d2phi_filtered_ar[i],
                     energy_den[i],pressure[i],r_e,coupling);

      source_mu[i] = source_mu_function(i, s_gp[i], r_e*mu_hat[i], r_e*nu_hat[i], r_e*phi_hat[i], dnu_filtered_ar[i], d2nu_filtered_ar[i],
                     dmu_filtered_ar[i], d2mu_filtered_ar[i], dphi_filtered_ar[i], d2phi_filtered_ar[i],
                     energy_den[i],pressure[i],r_e,coupling);

      source_phi[i] = source_phi_function(i, s_gp[i], dmu_filtered_ar[i],r_e);

      integral_mu1[i] = source_mu[i]/pow(1-s_gp[i],2);
      integral_mu2[i] = source_mu[i]*pow(1-s_gp[i],-1)/s_gp[i];
      integral_nu1[i] = source_nu[i]/pow(1-s_gp[i],2);
      integral_nu2[i] = source_nu[i]*pow(1-s_gp[i],-1)/s_gp[i];
      integral_phi[i] = source_phi[i]/(pow(1-s_gp[i],2)/sqrt(r_e));
    }

    for (i = 1; i < SDIV-1; ++i)
    {
      mu[i] = -(((1-s_gp[i])/s_gp[i])*trapz(slice(integral_mu1,0,i-1),slice(integral_mu1,0,i-1).size(), DS)
              + trapz(slice(integral_mu2,i,integral_mu2.size()-2), slice(integral_mu2,i,integral_mu2.size()-2).size(), DS));
      nu[i] = -(((1-s_gp[i])/s_gp[i])*trapz(slice(integral_nu1,0,i-1),slice(integral_nu1,0,i-1).size(), DS)
              + trapz(slice(integral_nu2,i,integral_mu2.size()-2), slice(integral_nu2,i,integral_mu2.size()-2).size(), DS));

      phi[i] = trapz(slice(integral_phi,0,i),slice(integral_phi,0,i).size(), DS);
    }

    mu[0] = -trapz(slice(integral_mu2,1,integral_mu2.size()-2), slice(integral_mu2,1,integral_mu2.size()-2).size(), DS);
    nu[0] = -trapz(slice(integral_nu2,1,integral_nu2.size()-2), slice(integral_nu2,1,integral_nu2.size()-2).size(), DS);
    phi[0] = trapz(slice(integral_phi,0,1),slice(integral_phi,0,1).size(), DS);

    mu[SDIV-1] = -((1-s_gp[SDIV-2])/s_gp[SDIV-2])*trapz(slice(integral_mu1,0,SDIV-2),slice(integral_mu1,0,SDIV-2).size(), DS);
    nu[SDIV-1] = -((1-s_gp[SDIV-2])/s_gp[SDIV-2])*trapz(slice(integral_nu1,0,SDIV-2),slice(integral_nu1,0,SDIV-2).size(), DS);
    phi[SDIV-1]= trapz(slice(integral_phi,0,SDIV-2),slice(integral_phi,0,SDIV-2).size(), DS);


    rel_error = (sqrt(r_e)-r_e_guess)/sqrt(r_e);

    iter+=1;

     if(print_option==0)
     {
       printf("---------------------------------\n");   
       std::cout<<"ITER ="<<" "<<iter<<'\n';
       std::cout<<"ERROR ="<<" "<<rel_error<<'\n';
       std::cout<<"r_e ="<<" "<<sqrt(r_e)<<'\n';
       printf("---------------------------------\n");   
     }

    if ((iter > 2)&&(fabs(rel_error)<= TOL || iter == MAX_ITER || isnan(r_e)))
    {
      break;
    }
    else{

      r_e_guess = sqrt(r_e);
      for (int m = 0; m < SDIV; ++m)
      {
        mu_guess[m] = RLX_FAC*mu[m] + (1-RLX_FAC)*mu_guess[m];
        nu_guess[m] = RLX_FAC*nu[m] + (1-RLX_FAC)*nu_guess[m];
        phi_guess[m] = RLX_FAC*phi[m] + (1-RLX_FAC)*phi_guess[m];
      }

    }
  }  
}


int main(int argc, char const *argv[])
{
	scanf("%s",eos_name);

	load_eos(eos_name);

  scanf("%lf",&coupling);
	scanf("%lf",&central_density);

  scanf("%lf",&TOL);
  scanf("%i",&MAX_ITER);
  scanf("%lf",&RLX_FAC);

  scanf("%i",&print_option);

  coupling/=KAPPA/pow(100000,2);

  make_center(central_density);
	make_grid();
	guess();
  main_iteration();

  std::vector<double> r_iso;

  for (auto s:s_gp)
    r_iso.push_back(sqrt(r_e)*s/(1-s));
  

  Mass = 2*r_iso[SDIV-3]*(sqrt(exp(mu[SDIV-3]))-1)*sqrt(KAPPA)*(C*C/G)/MSUN;
  Radius = sqrt(r_e)*exp(mu[int(SDIV/2)])*sqrt(KAPPA);

  switch(print_option){

      case 0:

      if(isnan(r_e)){
        std::cout<<"Failed to converge to a solution"<<'\n';
        break;
      }

      else{
        printf("\n");
        printf("  %s                  \n",eos_name);
        printf("  %6.5e  e_c             (10^15 gr/cm^3)\n",central_density);
        printf("\n");
        printf("  %.2e     alpha           (km^2)\n",coupling*KAPPA/pow(100000,2));
        printf("\n");
        printf("  %ix1        grid size              \n",SDIV);
        printf("  %2.2e     minimum relative error             \n",TOL);
        printf("  %2.2e     relative error             \n",fabs(rel_error));
        printf("  %d           maximum iterations             \n",MAX_ITER);
        printf("  %d            iteretions             \n",iter);
        printf("  %2.2e     relaxation factor             \n",RLX_FAC);
        printf("\n");
        printf("  %6.5e  M             (M_sun)\n",Mass);
        printf("  %6.5e  R             (km)\n",Radius/pow(10,5));
        printf("\n");
        break;

      }
      
      case 1:

     if(isnan(r_e)){
        std::cout<<"Failed to converge to a solution"<<'\n';
        break;
      }

      else{
        printf("\n");
        printf("  %s                  \n",eos_name);
        printf("  %6.5e  e_c             (10^15 gr/cm^3)\n",central_density);
        printf("\n");
        printf("  %.2e     alpha           (km^2)\n",coupling*KAPPA/pow(100000,2));
        printf("\n");
        printf("  %ix1        grid size              \n",SDIV);
        printf("  %2.2e     minimum relative error             \n",TOL);
        printf("  %2.2e     relative error             \n",fabs(rel_error));
        printf("  %d           maximum iterations             \n",MAX_ITER);
        printf("  %d            iteretions             \n",iter);
        printf("  %2.2e     relaxation factor             \n",RLX_FAC);
        printf("\n");
        printf("  %6.5e  M             (M_sun)\n",Mass);
        printf("  %6.5e  R             (km)\n",Radius/pow(10,5));
        printf("\n");
         
        printf("-------------------------------------------------------------------------------------\n");   
        printf("     s         r_is         mu           nu           phi       epsilon    pressure\n");
        printf("-------------------------------------------------------------------------------------\n");   


        for(size_t w =0;w<=SDIV;w++){

            printf("%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e \n",
             s_gp[w], r_iso[w], mu[w], nu[w], phi[w], energy_den[w], pressure[w]);
            }

        break;
    }

    return 0;
  }
}
