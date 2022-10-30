#include "../include/eos.hpp"
#include "../include/constants.hpp"
#include "../boost/math/interpolators/pchip.hpp"
#include <cmath>

using namespace SURFACE;
using namespace UNITS;
using namespace CGS;
using boost::math::interpolators::pchip; // Piecewise Cubic Hermite Interpolation

EOS::EOS(std::string eos_file)
{
  std::string eos;

  // set EoS path according to call type
  eos = EXECUTION::CALL_TYPE=="main_call" ? "EoS/" : "../EoS/";
  eos+=eos_file;

  // check to see if file is valid
  if((f_eos=fopen(eos.c_str(),"r")) == NULL ) {    
      std::cout<<"cannot open file "<<eos_file<<std::endl; 
      exit(0);
  }

  // load EoS file parameters in log vectors
  fscanf(f_eos,"%d\n",&n_tab);
  for(i=1;i<=n_tab;i++) {  
      fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ; 
      log_e_tab.push_back(log10(rho*C*C*KSCALE));
      log_p_tab.push_back(log10(p*KSCALE));
      log_h_tab.push_back(log10(h/(C*C)));
      log_n0_tab.push_back(log10(n0));
  }
  p_surface=pow(10,log_p_tab[0]);
  e_surface=pow(10,log_e_tab[0]);
  fclose(f_eos);
}

EOS::~EOS(){}

//----- pchip interpolation of tabulated EoS points -----//
double EOS::e_at_p(double pp) const
{
  auto spline = pchip<decltype(log_p_tab)>(log_p_tab,log_e_tab);
  double eTab = pp<p_surface? 0 : pow(10.0,spline(log10(pp)));
  return eTab;
}

double EOS::rho_at_p(double pp) const
{
  auto spline = pchip<decltype(log_p_tab)>(log_p_tab,log_rho_tab);
  double rhoTab = pp<p_surface ? 0.0 : pow(10.0,spline(log10(pp)));
  return rhoTab;
}

double EOS::p_at_e(double ee) const
{
    auto spline = pchip<decltype(log_e_tab)>(log_e_tab,log_p_tab);
    return pow(10.0,spline(log10(ee)));
}

double EOS::p_at_h(double hh) const
{ 
  auto spline = pchip<decltype(log_h_tab)>(log_h_tab,log_p_tab);
  double pTab = hh<=enthalpy_min ? 0.0 : pow(10.0,spline(log10(hh)));
  return pTab;
}

double EOS::h_at_p(double pp) const
{ 
    auto spline = pchip<decltype(log_p_tab)>(log_p_tab,log_h_tab);
    return pow(10.0,spline(log10(pp)));
}
//-----                                             -----//