#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;

    // Reionization parameters
    double zreion;
    double dzreion;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = Constants.x_start;
    const double x_end    = 0;
    
    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = 10000;
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;

    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // The two things we need to solve: Xe/ne and tau
    void solve_for_optical_depth_tau();

    // Solve the sound horizon ODE
    void solve_sound_horizon();

    // Splines contained in this class
    Spline log_Xe_of_x_spline{"Xe"};
    Spline log_Xe_reion_of_x_spline{"Xe_reion"};
    Spline tau_of_x_spline{"tau"}; 
    Spline tau_reion_of_x_spline{"tau_reion"};
    Spline g_tilde_of_x_spline{"g"};  
    Spline g_tilde_reion_of_x_spline{"g_reion"};
    Spline log_Xe_of_x_Saha_spline{"Xe_Saha"};
    Spline tau_of_x_Saha_spline{"tau_Saha"};
    Spline sound_horizon_spline{"Sound_Horizon"};

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp, double zreion, double dzreion);

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double tau_of_x(double x, bool reion) const;
    double dtaudx_of_x(double x, bool reion) const;
    double ddtauddx_of_x(double x, bool reion) const;
    double g_tilde_of_x(double x, bool reion) const;
    double dgdx_tilde_of_x(double x, bool reion) const;
    double ddgddx_tilde_of_x(double x, bool reion) const;
    double Xe_of_x(double x, bool reion) const;
    double ne_of_x(double x) const;
    double get_Yp() const;
    double get_zreion() const;
    double get_dzreion() const;
    double Xe_of_x_Saha(double x) const;
    double sound_horizon(double x) const;

    // Finding last scattering surface and recombination
    double last_scattering_surface() const;
    double last_scattering_surface_Saha() const;
    double x_at_recombination() const;
    double x_at_recombination_Saha() const;
};

#endif
