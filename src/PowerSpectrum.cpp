#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc,
    bool separate_components) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc),
  separate_components(separate_components)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  const double eta0 = cosmo -> eta_of_x(0);

  // Compute the spacings for the arrays such that the functions are properly sampled
  double delta_k_1 = (2.0*M_PI)/(16.0 * eta0);
  int npts_k_1 = std::trunc((k_max - k_min)/delta_k_1) + 1;
  Vector k_array_1 = Utils::linspace(k_min, k_max, npts_k_1);

  double delta_k_2 = (2.0*M_PI)/(6.0 * eta0);
  int npts_k_2 = std::trunc((k_max - k_min)/delta_k_2) + 1;
  Vector k_array_2 = Utils::linspace(k_min, k_max, npts_k_2);

  double delta_log_k = (2.0 * M_PI)/(6.0*eta0*k_max);
  int npts_log_k = std::trunc((log(k_max) - log(k_min))/delta_log_k) + 1;
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), npts_log_k);
  

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines(k_array_1);

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array_2);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //Get the Cell for each component of the temperature source function
  if(separate_components){
    auto cell_SW = solve_for_cell(log_k_array, SW_ell_of_k_spline, SW_ell_of_k_spline);
    cell_SW_spline.create(ells, cell_SW, "Cell_SW_of_ell");
    auto cell_ISW = solve_for_cell(log_k_array, ISW_ell_of_k_spline, ISW_ell_of_k_spline);
    cell_ISW_spline.create(ells, cell_ISW, "Cell_ISW_of_ell");
    auto cell_Doppler = solve_for_cell(log_k_array, Doppler_ell_of_k_spline, Doppler_ell_of_k_spline);
    cell_Doppler_spline.create(ells, cell_Doppler, "Cell_Doppler_of_ell");
    auto cell_Quad = solve_for_cell(log_k_array, Quad_ell_of_k_spline, Quad_ell_of_k_spline);
    cell_Quad_spline.create(ells, cell_Quad, "Cell_Quad_of_ell");
  }
  

  if(Constants.polarization){
    auto cell_EE = solve_for_cell(log_k_array, thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");
    auto cell_TE = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
    cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");
  }
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(Vector & k_array){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  
  double eta0 = cosmo -> eta_of_x(0);
  Vector k_eta_array;

  //Manually insert 0 for splining the Bessel function
  k_eta_array.push_back(0);

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    Vector bessel_array;

    //Manually insert 0 for splining the Bessel function
    bessel_array.push_back(0);
    auto append = [&] (double k){
      if(i == 0){
      k_eta_array.push_back(k * eta0);
      }

      bessel_array.push_back(Utils::j_ell(ell, k * eta0));
    };

    std::for_each(k_array.begin(), k_array.end(), append);

    j_ell_splines[i].create(k_eta_array, bessel_array, "j_ell");
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");
  
  
  double x_min = Constants.x_start;
  double x_max = 0;

  
  // Determine distance between points in the array to use trapezoid method of integration
  double delta_x = (x_max-x_min)/double(600-1);
  Vector x_array = Utils::linspace(x_min, x_max, 600);


  // Determine the particle horizon today
  double eta0 = cosmo -> eta_of_x(0);

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  for(size_t il = 0; il < ells.size(); il++){
    const int ell = ells[il];

    for(size_t ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // Initiate the integral result
      double integral = 0.0;
      
      for (size_t ix = 0; ix < x_array.size(); ix++){
        const double x = x_array[ix];
        double eta = cosmo -> eta_of_x(x);

        double bessel_factor = j_ell_splines[il](k*(eta0-eta));

        // Using the trapezoid method, multiply by 1/2 if dealing with either endpoint
        // (since they only enter the calculation once, unlike the remaining points)
        integral = (ix == 0 || ix == x_array.size()-1)? integral + delta_x/2.0 * bessel_factor * source_function(x, k) : integral + delta_x * bessel_factor * source_function(x, k);
      }

      result[il][ik] = integral;
    }
    
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //=======================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);
  for(size_t ell = 0; ell < thetaT_ell_of_k.size(); ell++){
    thetaT_ell_of_k_spline[ell].create(k_array, thetaT_ell_of_k[ell], "ThetaT_ell");
  }

  // Get individual contributions from each term in ST
  if(separate_components){
    SW_ell_of_k_spline = std::vector<Spline>(nells);
    ISW_ell_of_k_spline = std::vector<Spline>(nells);
    Doppler_ell_of_k_spline = std::vector<Spline>(nells);
    Quad_ell_of_k_spline = std::vector<Spline>(nells);
    std::function<double(double,double)> source_function_SW = [&](double x, double k){
      return pert->get_SW_term(x,k);
    };
    Vector2D SW_ell_of_k = line_of_sight_integration_single(k_array, source_function_SW);

    std::function<double(double,double)> source_function_ISW = [&](double x, double k){
      return pert->get_ISW_term(x,k);
    };
    Vector2D ISW_ell_of_k = line_of_sight_integration_single(k_array, source_function_ISW);

    std::function<double(double,double)> source_function_Doppler = [&](double x, double k){
      return pert->get_Doppler_term(x,k);
    };
    Vector2D Doppler_ell_of_k = line_of_sight_integration_single(k_array, source_function_Doppler);

    std::function<double(double,double)> source_function_Quad = [&](double x, double k){
      return pert->get_Quad_term(x,k);
    };
    Vector2D Quad_ell_of_k = line_of_sight_integration_single(k_array, source_function_Quad);

    for(size_t ell = 0; ell < thetaT_ell_of_k.size(); ell++){
      SW_ell_of_k_spline[ell].create(k_array, SW_ell_of_k[ell], "SW_ell");
      ISW_ell_of_k_spline[ell].create(k_array, ISW_ell_of_k[ell], "ISW_ell");
      Doppler_ell_of_k_spline[ell].create(k_array, Doppler_ell_of_k[ell], "Doppler_ell");
      Quad_ell_of_k_spline[ell].create(k_array, Quad_ell_of_k[ell], "Quad_ell");
    }
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){
    thetaE_ell_of_k_spline = std::vector<Spline>(nells);

    std::function<double(double,double)> source_function_E = [&](double x, double k){
      return pert->get_Source_E(x,k);
    };

    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, source_function_E);

    for(size_t ell = 0; ell < thetaE_ell_of_k.size(); ell++){
      thetaE_ell_of_k_spline[ell].create(k_array, thetaE_ell_of_k[ell], "ThetaT_ell");
    }
  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
  
  Vector k_array = exp(log_k_array); 
  double delta_log_k = (log_k_array.back() - log_k_array[0])/(log_k_array.size());
  Vector result(nells, 0.0);

  for(size_t ell = 0; ell < result.size(); ell++){

    for(size_t ik = 1; ik < log_k_array.size(); ik++){

      result[ell] = result[ell] + 4 * M_PI * (delta_log_k/2) 
        * (f_ell_spline[ell](k_array[ik])*g_ell_spline[ell](k_array[ik])*primordial_power_spectrum(k_array[ik])
        + f_ell_spline[ell](k_array[ik-1])*g_ell_spline[ell](k_array[ik-1])*primordial_power_spectrum(k_array[ik-1]));
    }

  }

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  double pofk = 0.0;

  // Quantities needed to compute the gauge invariant matter density perturbation (a is not needed since it is 1 for x=0)
  const double c = Constants.c;
  const double H0 = cosmo -> get_H0();
  double OmegaM = cosmo->get_OmegaB(0) + cosmo->get_OmegaCDM(0);
  double Phi = pert -> get_Phi(0, k);

  // Gauge invariant matter density perturbation 
  double Delta_M = (pow(c*k, 2.0) * Phi)/(1.5 * OmegaM * pow(H0, 2.0));

  pofk = pow(Delta_M, 2.0)*(2.0 * pow(M_PI, 2.0))/pow(k, 3.0) * primordial_power_spectrum(k);

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}
double PowerSpectrum::get_k_eq() const{
  // Fetch the relevant quantities (x,a and H) at matter radiation equality
  const double x_eq = cosmo -> matter_radiation_equality();
  const double a_eq = exp(x_eq);
  const double H_eq = cosmo -> H_of_x(x_eq);

  return a_eq*H_eq/Constants.c;
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename1, std::string filename2) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename1.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);

  // Start printing matter power spectrum and transfer functions

  std::ofstream fp2(filename2.c_str());
  double eta0 = cosmo -> eta_of_x(0);
  auto print_data_k = [&] (const double k) {
    fp2 << k << " ";
    fp2 << get_matter_power_spectrum(0,k)<< " ";
    for (size_t ell = 0; ell<ells.size(); ell++){
      fp2 << thetaT_ell_of_k_spline[ell](k) << " ";
    }
    fp2 << "\n";
  };
  double delta_k = (2.0*M_PI)/(6.0 * eta0);
  int npts_k = std::trunc((k_max - k_min)/delta_k) + 1;
  Vector k_array = Utils::linspace(k_min, k_max, npts_k);
  std::for_each(k_array.begin(), k_array.end(), print_data_k);
  fp2 << "k at matter radiation equality:" << "\n";
  fp2 << get_k_eq();

  if(separate_components){
    std::string filename3 = "contributions.txt";
    std::ofstream fp3(filename3.c_str());
    auto print_data_cpts = [&] (const double ell) {
      double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
      fp3 << ell                                 << " ";
      fp3 << cell_SW_spline( ell ) * normfactor  << " ";
      fp3 << cell_ISW_spline( ell ) * normfactor  << " ";
      fp3 << cell_Doppler_spline( ell ) * normfactor  << " ";
      fp3 << cell_Quad_spline( ell ) * normfactor  << " ";
      fp3 << "\n";
    };
    std::for_each(ellvalues.begin(), ellvalues.end(), print_data_cpts);
  }
}

