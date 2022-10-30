#ifndef KEH_CST_HPP
#define KEH_CST_HPP

#include "eos.hpp"
#include "constants.hpp"
#include <vector>

using String = std::string;
using Vector = std::vector<double>;
using namespace INTEGRATION;

// This is the main class that incorporates the Komatsu-Eriguchi-Hachisu/Cook-Shapiro-Teukolsky method // 
class KEH_CST
{
    public:

        KEH_CST(String eosInput, double couplingInput, double pressureInput, double accuracyInput,
                double relaxationInput, int maximumIterationInput, int printOptionInput);

    private:

        void make_grid(void); // create computational grid
        void make_center(double); // computes values at the starting point of the grid
        void load_guess(String eosInput, double couplingInput, double pressureInput); // loads guess solution files      
        Vector slice_vec(Vector const&, int, int) const; //vector slicer used in main routine integration
        Vector linspace(double start_in, double end_in, int num_in) const; // linearly spaced vector
        double firstOrderDerivative(int, double*, double) const; // 
        double trapz(Vector, double, double) const; // trapezoidal rule
        void relaxationFactor(double&); // relaxation factor adaptor
        void main_iteration(String eosInput, double couplingInput, double pressureInput); // main iterative routine
        
        /* source functions */
        double source_nu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				                          double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				                          double e_aw, double p_aw, double r_e, double a) const;
        double source_mu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
			                          	double dnds2, double dmds, double dmds2, double dfds, double dfds2,
			                          	double e_aw, double p_aw, double r_e, double a) const;
        double source_phi_function(int i, double s, double dmds, double r_e) const;
        /*                   */
    
    private:

        double central_pressure, central_density, central_enthalpy;
        double s_gp[SDIV]; // grid points
        double pressure[SDIV], energy_den[SDIV], enthalpy[SDIV]; // pressure, energy density and enthalpy points
        double source_mu[SDIV], source_nu[SDIV], source_phi[SDIV]; // source points
        double mu[SDIV],nu[SDIV],phi[SDIV]; // metric function and scalar points
        double mu_hat[SDIV], mu_guess[SDIV], nu_hat[SDIV], nu_guess[SDIV], phi_hat[SDIV], phi_guess[SDIV]; // re-factored and guess, metric and scalar points
        Vector dmu,dnu,dphi,d2mu,d2nu,d2phi; // derivative vectors
        Vector dmu_filtered, dnu_filtered, dphi_filtered, d2mu_filtered, d2nu_filtered, d2phi_filtered; // Savitzky-Golay filtered derivative vectors
        Vector integral_mu1, integral_mu2,integral_nu1,integral_nu2,integral_phi; // integrals used in main integration routine
        Vector r_iso; // isotropic coordinate vector
        double nu_hat_eq; // re-factored nu metric at surface
        double r_e, r_e_guess; // coordinate radius at surface and its guess
        double rel_error; // relative error in concurrent r_e evaluations
        double coupling; // coupling value
        double relaxation, accuracy; // relaxation factor and accuracy in r_e 
        int iter, max_iter; // iteration and maximum iterations
        int convergence_check; // checks if method converges
        int print_option;
        std::string eos_name;
        EOS eos;
    
    public:

        void compute_MR(); // computes model
        void printModel() const; // prints model 
        double mass,radius,relaxation_output;
        int convergence_check_output;
};
#endif