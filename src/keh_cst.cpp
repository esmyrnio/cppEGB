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

KEH_CST::KEH_CST(eosInput eos_name, couplingInput coupling, pressureInput central_pressure,
                 accuracyInput accuracy, relaxationInput relaxation, maximumIterationInput max_iter, printOptionInput print_option)
                 :eos(eos_name), convergence_check(1), dmu(SDIV), dnu(SDIV), dphi(SDIV), d2mu(SDIV), d2nu(SDIV), d2phi(SDIV), integral_mu1(SDIV), integral_mu2(SDIV),
                 integral_nu1(SDIV), integral_nu2(SDIV), integral_phi(SDIV)
                {
                    this->eos_name = eos_name.c_str();
                    this->coupling = coupling*pow(100000,2)/KAPPA;
                    this->central_pressure = central_pressure*pow(10,35)*KSCALE;
                    this->accuracy = accuracy;
                    this->relaxation = relaxation;
                    this->max_iter = max_iter;
                    this->print_option = print_option;

                    central_density = eos.e_at_p(this->central_pressure);
                    central_enthalpy = eos.h_at_p(this->central_pressure);

                    for(size_t s=1;s<=SDIV;s++) {
                        s_gp[s-1]=(1.0*s-1.0)/(SDIV-1);
                    }
                 }

double KEH_CST::firstOrderDerivative(int i, double* f, double crit){
  if (i==0) {return crit*(f[i+1]-f[i])/DS;}
  else if (i==SDIV-1)  {return crit*(f[i]-f[i-1])/DS;}  
  else {return crit*(f[i+1]-f[i-1])/(2*DS);}
}

double KEH_CST::trapz(std::vector<double> ar, double size, double step){
  double res;
  double sum = 0;
  for(double x:ar)
    sum+=x;
  res = step*(sum+(ar[0]+ar[size-1])/2);
  return res;
}

std::vector<double> KEH_CST::slice_vec(std::vector<double> const& v, int X, int Y)
{
    auto first = v.begin() + X;
    auto last = v.begin() + Y + 1;
    std::vector<double> vector(first, last);
    return vector;
}

std::vector<double> KEH_CST::linspace(double start_in, double end_in, int num_in)
{
  std::vector<double> linspaced;
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

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

void KEH_CST::load_guess(eosInput eos_input, couplingInput coupling_input, pressureInput pressure_input){

  std::string eos_dir;
  std::regex target(".rns1.1.txt");
  eos_dir = std::regex_replace(eos_input,target,"");

  std::vector<double> pressures = linspace(0.1,30.0,30);
  std::vector<double> couplings = linspace(1.0,70.0,10);

  std::for_each(std::begin(pressures),std::end(pressures), [&pressure_input](double& d_p) { d_p=pow(d_p-pressure_input,2); });
  std::for_each(std::begin(couplings),std::end(couplings), [&coupling_input](double& d_a) { d_a=pow(d_a-coupling_input,2); });

  auto it_pressure = std::min_element(std::begin(pressures),std::end(pressures));
  auto it_coupling = std::min_element(std::begin(couplings),std::end(couplings));

  int p_nearest = std::distance(std::begin(pressures),it_pressure); 
  int a_nearest = std::distance(std::begin(couplings),it_coupling);

  FILE *guess_file;

  std::string guess_file_name;
    
  if (call_type=="get_MR")
  {
    guess_file_name="../guess_solutions/";
  }
  else if (call_type=="main_call")
  {
    guess_file_name="guess_solutions/";
  }

  std::stringstream p_extension;
  std::stringstream a_extension;

  std::vector<double> pressures_cp = linspace(0.1,30.0,30);
  std::vector<double> couplings_cp = linspace(1.0,70.0,10);
  a_extension<<"alpha_"<<std::setprecision(4)<<couplings_cp[a_nearest];
  p_extension<<"pressure_"<<std::setprecision(4)<<pressures_cp[p_nearest];
  
  guess_file_name = guess_file_name+eos_dir+'/'+a_extension.str()+'/'+p_extension.str();

  if((guess_file=fopen(guess_file_name.c_str(),"r")) == NULL ) { 
      std::cout<<"Cannot open file "<<guess_file_name<<std::endl;   
      exit(0);
  }

  double mu,nu,phi;
  fscanf(guess_file,"%lf\n",&r_e_guess);
  for(int i=0;i<SDIV;i++){
    fscanf(guess_file,"%lf %lf %lf\n",&mu,&nu,&phi);
    mu_guess[i]=mu;
    nu_guess[i]=nu;
    phi_guess[i]=phi;
  }

}

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

double KEH_CST::source_phi_function(int i, double s, double dmds, double r_e){
	double ss = pow(1-s,2)/sqrt(r_e);
	double dm = ss*dmds;
	return dm;
}

void KEH_CST::relaxationFactor(double& _relaxationFactor){
  _relaxationFactor-=0.01;
  if(_relaxationFactor<1e-2){
  convergence_check=0;
  _relaxationFactor=1e-2;
  }
}

void KEH_CST::main_iteration(eosInput _eos_name, couplingInput _coupling, pressureInput _pressure){

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
      if(abs(enthalpy[i])-abs(enthalpy[i-1])>1e-4) break;
      pressure[i] = eos.p_at_h(enthalpy[i]);
      energy_den[i] = eos.e_at_p(pressure[i]);

      dmu[i] = firstOrderDerivative(i,mu_hat,r_e);
      dnu[i] = firstOrderDerivative(i,nu_hat,r_e);
      dphi[i] = firstOrderDerivative(i,phi_hat,r_e);
    }
 
    if((i<SDIV-1) && (convergence_check!=0)){
      relaxationFactor(relaxation);
      load_guess(_eos_name,_coupling,_pressure);
      main_iteration(_eos_name,_coupling,_pressure);
    }

    dmu_filtered = sg_smooth(dmu,31,3);

    for (i = 0; i < SDIV; ++i)
    {
      d2mu[i] = firstOrderDerivative(i,&dmu_filtered[0],1.0);
      d2nu[i] = firstOrderDerivative(i,&dnu[0],1.0);
      d2phi[i] = firstOrderDerivative(i,&dphi[0],1.0);
    }

    d2mu_filtered = sg_smooth(d2mu,31,3);

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

    rel_error = abs((r_e_guess-sqrt(r_e))/sqrt(r_e));

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
