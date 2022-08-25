#ifndef KEH_CST_HPP
#define KEH_CST_HPP

#include "eos.hpp"
#include "constants.hpp"
#include <vector>

typedef std::string eosInput;
typedef double couplingInput;
typedef double pressureInput;
typedef double toleranceInput;
typedef double accuracyInput;
typedef double relaxationInput;
typedef int maximumIterationInput;
typedef int printOptionInput;

using namespace INTEGRATION;

class KEH_CST{

    public:

        KEH_CST(eosInput, couplingInput, pressureInput, accuracyInput, relaxationInput,maximumIterationInput,printOptionInput);

    private:
        void make_center(double);
        void make_grid(void);
        void load_guess(eosInput, couplingInput, pressureInput);
        
        std::vector<double> slice_vec(std::vector<double> const&, int, int);
        std::vector<double> linspace(double start_in, double end_in, int num_in);

        static double firstOrderDerivative(int, double*, double);
        static double trapz(std::vector<double>, double, double);

        static double source_nu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				double e_aw, double p_aw, double r_e, double a);
        static double source_mu_function(int i, double s, double mu_t, double nu_t, double phi_t, double dnds,
				double dnds2, double dmds, double dmds2, double dfds, double dfds2,
				double e_aw, double p_aw, double r_e, double a);
        static double source_phi_function(int i, double s, double dmds, double r_e);

        void relaxationFactor(double&);
        void main_iteration(eosInput, couplingInput, pressureInput);

    private:
        double central_pressure, central_density, central_enthalpy;
        double s_gp[SDIV];
        double pressure[SDIV], energy_den[SDIV], enthalpy[SDIV];
        double source_mu[SDIV], source_nu[SDIV], source_phi[SDIV];
        double mu[SDIV],nu[SDIV],phi[SDIV];
        double mu_hat[SDIV], mu_guess[SDIV], nu_hat[SDIV], nu_guess[SDIV], phi_hat[SDIV], phi_guess[SDIV];
        std::vector<double> dmu,dnu,dphi,
                            d2mu,d2nu,d2phi;
        std::vector<double> dmu_filtered, dnu_filtered, dphi_filtered, d2mu_filtered, 
                            d2nu_filtered, d2phi_filtered;
        std::vector<double> integral_mu1, integral_mu2,integral_nu1,integral_nu2,integral_phi;
        std::vector<double> r_iso;
        double nu_hat_eq;
        double r_e, r_e_guess;
        double rel_error;
        double coupling, relaxation, accuracy;
        int iter, max_iter;
        int convergence_check;
        int print_option;
        std::string eos_name;
        EOS eos;

    public:
        void compute_MR();
        void printModel();
        double mass,radius,relaxation_output;
        int convergence_check_output;



};

#endif
