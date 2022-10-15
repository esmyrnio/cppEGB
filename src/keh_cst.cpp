#include "../include/keh_cst.hpp"
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../include/savgol.hpp"
#include <fstream> 
#include <math.h>

using namespace EXECUTION;
using namespace CGS;
using namespace UNITS;

KEH_CST::KEH_CST(String eos_name, double coupling, double central_pressure,
                 double accuracy, double relaxation, int max_iter, int print_option)
                 :eos_name(eos_name.c_str()),coupling(coupling*pow(100000,2)/KAPPA),central_pressure(central_pressure*pow(10,35)*KSCALE),
                  accuracy(accuracy),relaxation(relaxation),max_iter(max_iter),print_option(print_option),
                  eos(eos_name), convergence_check(1), dmu(SDIV), dnu(SDIV), dphi(SDIV), d2mu(SDIV), d2nu(SDIV), d2phi(SDIV), integral_mu1(SDIV), integral_mu2(SDIV),
                  integral_nu1(SDIV), integral_nu2(SDIV), integral_phi(SDIV)
                {
                    central_density = eos.e_at_p(this->central_pressure);
                    central_enthalpy = eos.h_at_p(this->central_pressure);
                    // construct grid
                    for(size_t s=1;s<=SDIV;s++) {
                        s_gp[s-1]=(1.0*s-1.0)/(SDIV-1);
                    }
                 }
// returns first order derivate
double KEH_CST::firstOrderDerivative(int i, double* f, double crit){
  if (i==0) {return crit*(f[i+1]-f[i])/DS;}
  else if (i==SDIV-1)  {return crit*(f[i]-f[i-1])/DS;}  
  else {return crit*(f[i+1]-f[i-1])/(2*DS);}
}
// simple trapezoidal integration
double KEH_CST::trapz(Vector ar, double size, double step){
  double res;
  double sum = 0;
  for(double x:ar)
    sum+=x;
  res = step*(sum+(ar[0]+ar[size-1])/2);
  return res;
}
// slices vector at given points
Vector KEH_CST::slice_vec(Vector const& v, int X, int Y)
{
    auto first = v.begin() + X;
    auto last = v.begin() + Y + 1;
    Vector vector(first, last);
    return vector;
}
// returns a linearly spaced vector from "start" to "end" with "num" number of points
Vector KEH_CST::linspace(double start, double end, int num)
{
  Vector linspaced;
  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }
  double delta = (end - start) / (num - 1);
  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); 
  return linspaced;
}
/* loads the GR solutions to be used as guess in the main method. This part is used for optimization of the KEH/CST method and performance.
 Alternatively one can simply integrate the TOV system in isotropic coordinates and used the output as guess */
void KEH_CST::load_guess(String eos_input, double coupling_input, double pressure_input){

  /* the guess solutions are stored w.r.t the equation of state, coupling (ranging from 1.0 to 70.0 km^2)
     and central pressure (ranging from 0.1 to 30.0 dyn/cm^2). So in order to load them we must follow the same format. */
  String eos_dir;
  std::regex target(".rns1.1.txt");
  eos_dir = std::regex_replace(eos_input,target,"");
  
  /* grid used when calculating the guess solutions */
  Vector pressures = linspace(0.1,30.0,30);
  Vector couplings = linspace(1.0,70.0,10);

  /* find the closest point of the grid to input pressure and coupling */
  std::for_each(std::begin(pressures),std::end(pressures), [&pressure_input](double& d_p) { d_p=pow(d_p-pressure_input,2); });
  std::for_each(std::begin(couplings),std::end(couplings), [&coupling_input](double& d_a) { d_a=pow(d_a-coupling_input,2); });
  auto it_pressure = std::min_element(std::begin(pressures),std::end(pressures));
  auto it_coupling = std::min_element(std::begin(couplings),std::end(couplings));
  int p_nearest = std::distance(std::begin(pressures),it_pressure); 
  int a_nearest = std::distance(std::begin(couplings),it_coupling);

  /* construct the appropriate extension for coupling and pressure */
  Vector pressures_cp = linspace(0.1,30.0,30);
  Vector couplings_cp = linspace(1.0,70.0,10);
  std::stringstream p_extension,a_extension;
  a_extension<<"alpha_"<<std::setprecision(4)<<couplings_cp[a_nearest];
  p_extension<<"pressure_"<<std::setprecision(4)<<pressures_cp[p_nearest];
  
  /* load the guess solution file */
  FILE *guess_file;
  String guess_file_name;
  guess_file_name = CALL_TYPE=="get_MR" ? "../guess_solutions/" : "guess_solutions/";  
  guess_file_name = guess_file_name+eos_dir+'/'+a_extension.str()+'/'+p_extension.str();
  if((guess_file=fopen(guess_file_name.c_str(),"r")) == NULL ) { 
      std::cout<<"Cannot open file "<<guess_file_name<<std::endl;   
      exit(0);
  }
  double mu,nu,phi;
  // store the file's guess solution in guess solution arrays
  fscanf(guess_file,"%lf\n",&r_e_guess);
  for(int i=0;i<SDIV;i++){
    fscanf(guess_file,"%lf %lf %lf\n",&mu,&nu,&phi);
    mu_guess[i]=mu;
    nu_guess[i]=nu;
    phi_guess[i]=phi;
  }
}
// returns the source for the nu metric function at given point 
double KEH_CST::source_nu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				double e_aw, double p_aw, double r_e, double a){

                double c = 4*M_PI;
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
                if(i!=0)

                {return r*(emu*(c*emu*r2*(3*p+eps)-r*dn*(2+r*(dm+dn))+a*memu*
                            (dphi*(2*r2*pow(dm,3)-2*r2*pow(dm,2)*(dn+dphi)-2*r*pow(dn,2)
                            *(-2+r*dphi)-r*(pow(dphi,2)*(-4+r*dphi)+4*ddm)+2*dm
                            *(-2+r*(2*r*pow(dn,2)+dn*(4-r*dphi)+dphi*(-2+r*dphi)-2*r*ddm))
                            +2*dn*(2+r*(dphi*(-4+r*dphi)+2*r*ddm)))+2*r*(-r*pow(dm,2)+2*dm*
                            (-1+r*dn)+r*pow(dphi,2)+dn*(2-2*r*dphi))*ddphi)))
                            /(emu*r+2*a*dphi*(-2-2*r*dm+r*dphi)) + term ;}
                    
                else return 0;
}
// returns the source for the mu metric function at given point
double KEH_CST::source_mu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				double e_aw, double p_aw, double r_e, double a){

                    double c = 8*M_PI;
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

                    if(i!=0)
                    
                        {return -r*((c*pow(emu,2)*r2*eps+emu*r*dm*(4+r*dm)+a*(2*dm-dphi)
                        *dphi*(-4+r2*(2*pow(dm,2)-2*dm*dphi+pow(dphi,2)))
                        -4*a*r*(dm-dphi)*(2+r*dm-r*dphi)*ddphi)/
                        (2*(emu*r+2*a*dphi*(-2-2*r*dm+r*dphi)))) + term ;}
                    
                    else return 0; 
}
// returns the source for the scalar at given point 
double KEH_CST::source_phi_function(int i, double s, double dmds, double r_e){
	double ss = pow(1-s,2)/sqrt(r_e);
	double dm = ss*dmds;
	return dm;
}
// routine to adapt the method's relaxation factor for convergence, whenever needed
void KEH_CST::relaxationFactor(double& _relaxationFactor){
  _relaxationFactor-=0.01;
  if(_relaxationFactor<1e-2){
  convergence_check=0;
  _relaxationFactor=1e-2;
  }
}
// main routine
void KEH_CST::main_iteration(String _eos_name, double _coupling, double _pressure){
  iter=0;
  std::size_t i;
  while(true){
    for (i = 0; i < SDIV; ++i)
    {
      mu_hat[i] = mu_guess[i]/pow(r_e_guess,2);
      nu_hat[i] = nu_guess[i]/pow(r_e_guess,2);
      phi_hat[i] = phi_guess[i]/pow(r_e_guess,2);
    }
    nu_hat_eq = nu_hat[int(SDIV/2)];
    r_e = central_enthalpy/(nu_hat_eq-nu_hat[0]);
    pressure[0]=central_pressure;
    energy_den[0]=central_density;
    enthalpy[0]=central_enthalpy;
    for (i = 1; i < SDIV; ++i)
    {
      enthalpy[i]= SURFACE::enthalpy_min + r_e*(nu_hat_eq-nu_hat[i]);
      if(i>int(SDIV/2)) enthalpy[i]=0.0;
      if(abs(enthalpy[i])-abs(enthalpy[i-1])>1e-4) break; // checks if enthalpy starts increasing, in which case the method will not converge to the right solution
      pressure[i] = eos.p_at_h(enthalpy[i]);
      energy_den[i] = eos.e_at_p(pressure[i]);
      dmu[i] = firstOrderDerivative(i,mu_hat,r_e);
      dnu[i] = firstOrderDerivative(i,nu_hat,r_e);
      dphi[i] = firstOrderDerivative(i,phi_hat,r_e);
    }
    // convergence check
    if((i<SDIV-1) && (convergence_check!=0)){
      relaxationFactor(relaxation);
      load_guess(_eos_name,_coupling,_pressure);
      main_iteration(_eos_name,_coupling,_pressure);
    }
    dmu_filtered = sg_smooth(dmu,31,3); // filter the first order derivative of the mu metric function to remove oscillations
    for (i = 0; i < SDIV; ++i)
    {
      d2mu[i] = firstOrderDerivative(i,&dmu_filtered[0],1.0);
      d2nu[i] = firstOrderDerivative(i,&dnu[0],1.0);
      d2phi[i] = firstOrderDerivative(i,&dphi[0],1.0);
    }
    d2mu_filtered = sg_smooth(d2mu,31,3); // filter the first order derivative of the mu metric function to remove oscillations
    for (i = 0; i < SDIV; ++i)
    {
      source_nu[i] = source_nu_function(i, s_gp[i], r_e*mu_hat[i], r_e*nu_hat[i], r_e*phi_hat[i], dnu[i], d2nu[i],
                     dmu_filtered[i], d2mu_filtered[i], dphi[i], d2phi[i],
                     energy_den[i],pressure[i],r_e,coupling);
      source_mu[i] = source_mu_function(i, s_gp[i], r_e*mu_hat[i], r_e*nu_hat[i], r_e*phi_hat[i], dnu[i], d2nu[i],
                     dmu_filtered[i], d2mu_filtered[i], dphi[i], d2phi[i],
                     energy_den[i],pressure[i],r_e,coupling);
      source_phi[i] = source_phi_function(i, s_gp[i], dmu_filtered[i],r_e);
      integral_mu1[i] = source_mu[i]/pow(1-s_gp[i],2);
      integral_mu2[i] = source_mu[i]*pow(1-s_gp[i],-1)/s_gp[i];
      integral_nu1[i] = source_nu[i]/pow(1-s_gp[i],2);
      integral_nu2[i] = source_nu[i]*pow(1-s_gp[i],-1)/s_gp[i];
      integral_phi[i] = source_phi[i]/(pow(1-s_gp[i],2)/sqrt(r_e));
    }
    for (i = 1; i < SDIV-1; ++i)
    {
      mu[i] = -(((1-s_gp[i])/s_gp[i])*trapz(slice_vec(integral_mu1,0,i-1),slice_vec(integral_mu1,0,i-1).size(), DS)
              + trapz(slice_vec(integral_mu2,i,integral_mu2.size()-2), slice_vec(integral_mu2,i,integral_mu2.size()-2).size(), DS));
      nu[i] = -(((1-s_gp[i])/s_gp[i])*trapz(slice_vec(integral_nu1,0,i-1),slice_vec(integral_nu1,0,i-1).size(), DS)
              + trapz(slice_vec(integral_nu2,i,integral_mu2.size()-2), slice_vec(integral_nu2,i,integral_mu2.size()-2).size(), DS));
      phi[i] = trapz(slice_vec(integral_phi,0,i),slice_vec(integral_phi,0,i).size(), DS);
    }
    mu[0] = -trapz(slice_vec(integral_mu2,1,integral_mu2.size()-2), slice_vec(integral_mu2,1,integral_mu2.size()-2).size(), DS);
    nu[0] = -trapz(slice_vec(integral_nu2,1,integral_nu2.size()-2), slice_vec(integral_nu2,1,integral_nu2.size()-2).size(), DS);
    phi[0] = trapz(slice_vec(integral_phi,0,1),slice_vec(integral_phi,0,1).size(), DS);
    mu[SDIV-1] = -((1-s_gp[SDIV-2])/s_gp[SDIV-2])*trapz(slice_vec(integral_mu1,0,SDIV-2),slice_vec(integral_mu1,0,SDIV-2).size(), DS);
    nu[SDIV-1] = -((1-s_gp[SDIV-2])/s_gp[SDIV-2])*trapz(slice_vec(integral_nu1,0,SDIV-2),slice_vec(integral_nu1,0,SDIV-2).size(), DS);
    phi[SDIV-1]= trapz(slice_vec(integral_phi,0,SDIV-2),slice_vec(integral_phi,0,SDIV-2).size(), DS);
    rel_error = abs((r_e_guess-sqrt(r_e))/sqrt(r_e)); // relative error that controls iterations
    if ((iter >= 2)&&(rel_error<= accuracy) || convergence_check==0)
    { break;
    }
    else if(isnan(r_e) || (iter==max_iter)&&(rel_error>10*accuracy)&&(convergence_check!=0)){
      relaxationFactor(relaxation);
      load_guess(_eos_name,_coupling,_pressure);
      main_iteration(_eos_name,_coupling,_pressure);
      }
    else{
      iter+=1;
      r_e_guess = sqrt(r_e);
      for (int m = 0; m < SDIV; ++m)
      {
        mu_guess[m] = relaxation*mu[m] + (1-relaxation)*mu_guess[m];
        nu_guess[m] = relaxation*nu[m] + (1-relaxation)*nu_guess[m];
        phi_guess[m] = relaxation*phi[m] + (1-relaxation)*phi_guess[m];
      }
    }
  }  
}
// computes model's gravitational mass and radius
void KEH_CST::compute_MR(){
  load_guess(eos_name,coupling*KAPPA/pow(100000,2), central_pressure/pow(10,35)/KSCALE);
  main_iteration(eos_name,coupling*KAPPA/pow(100000,2), central_pressure/pow(10,35)/KSCALE);
  for (auto s:s_gp)
    r_iso.push_back(sqrt(r_e)*s/(1-s));
  mass = 2*r_iso[SDIV-2]*(sqrt(exp(mu[SDIV-2]))-1)*sqrt(KAPPA)*(C*C/G)/MSUN;
  radius = sqrt(r_e)*exp(mu[int(SDIV/2)])*sqrt(KAPPA)/pow(10,5);
  relaxation_output=relaxation;
  convergence_check_output=convergence_check;
}
void KEH_CST::printModel(){
  switch (print_option)
  {
  case 0:
      if(std::isnan(r_e)){
          std::cout<<"Failed to converge to a solution"<<'\n';
        }
        else{
          printf("\n");
          printf("  %s                  \n",eos_name.c_str());
          printf("  %6.5e  p_c             (10^35 dyn/cm^2)\n",central_pressure/pow(10,35)/KSCALE);
          printf("\n");
          printf("  %lf     alpha           (km^2)\n",coupling*KAPPA/pow(100000,2));
          printf("\n");
          printf("  %ix1        grid size              \n",SDIV);
          printf("  %2.2e     minimum relative error             \n",accuracy);
          printf("  %2.2e     relative error             \n",abs(rel_error));
          printf("  %d           maximum iterations             \n",max_iter);
          printf("  %d            iteretions             \n",iter);
          printf("  %2.2e     relaxation factor             \n",relaxation);
          printf("\n");
          printf("  %6.5e  M             (M_sun)\n",mass);
          printf("  %6.5e  R             (km)\n",radius);
          printf("\n");

        }
    break;
  case 1:
    if(std::isnan(r_e)){
          std::cout<<"Failed to converge to a solution"<<'\n';
        }
    else{
      printf("\n");
      printf("  %s                  \n",eos_name.c_str());
      printf("  %6.5e  p_c             (10^35 dyn/cm^2)\n",central_pressure/pow(10,35)/KSCALE);
      printf("\n");
      printf("  %lf     alpha           (km^2)\n",coupling*KAPPA/pow(100000,2));
      printf("\n");
      printf("  %ix1        grid size              \n",SDIV);
      printf("  %2.2e     minimum relative error             \n",accuracy);
      printf("  %2.2e     relative error             \n",abs(rel_error));
      printf("  %d           maximum iterations             \n",max_iter);
      printf("  %d            iteretions             \n",iter);
      printf("  %2.2e     relaxation factor             \n",relaxation);
      printf("\n");
      printf("  %6.5e  M             (M_sun)\n",mass);
      printf("  %6.5e  R             (km)\n",radius);
      printf("\n");

    }
    for(int i=0;i<SDIV;i++){
      printf("%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e \n",
        s_gp[i], r_iso[i], mu[i], nu[i], phi[i], energy_den[i], pressure[i]);
    }
  }
}
