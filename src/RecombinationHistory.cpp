#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp, double zreion, double dzreion) :
  cosmo(cosmo),
  Yp(Yp),
  zreion(zreion),
  dzreion(dzreion)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  // Compute and spline the sound horizon
  solve_sound_horizon();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr;
  Vector ne_arr;
  Vector Xe_reion_arr;
  Vector Xe_Saha_arr;

  // Used to save the index of x_array at which the Saha regime ends
  int splitting_index;

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime? Only check the splitting index if we were
    // in the Saha regime to begin with
    if(saha_regime && Xe_current < Xe_saha_limit){
      saha_regime = false;
      splitting_index = i;
    };

    if(saha_regime){
     Xe_arr.push_back(log(Xe_current));
     Xe_reion_arr.push_back(log(Xe_current));
     ne_arr.push_back(ne_current); 
     Xe_Saha_arr.push_back(log(Xe_current));
    } else { 
      // Even if no longer in Saha regime, add points to the array that
      // computes the Saha prediction for all x to compare with peebles equation
     if (log(Xe_current) == log(Xe_current)){ 
        Xe_Saha_arr.push_back(log(Xe_current));
     } else {
        // If the log of the Saha approximation returns a NaN due to an underflow to 0, simply fill the
        // rest of the array with the last valid value (by all intents and purposes, that is 0)
        for(int j = i; j < npts_rec_arrays; j++){
          Xe_Saha_arr.push_back(Xe_Saha_arr.back());
        };
        break;
     };
      };
  }
  ODESolver peebles_Xe_ode;
  ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
    return rhs_peebles_ode(x, Xe, dXedx);
  };

  double Xe0 = exp(Xe_arr.back());
  Vector Xe_i = {Xe0};
  Vector new_x_array;
  for(size_t i = splitting_index; i < x_array.size(); i++){
    new_x_array.push_back(x_array[i]);
  };

  peebles_Xe_ode.solve(dXedx, new_x_array, Xe_i, gsl_odeiv2_step_rk4);

  auto Xe_arr_peebles = peebles_Xe_ode.get_data_by_component(0);

  for(size_t i = 0; i < Xe_arr_peebles.size(); i++){
    // Compute the correction provided by reionization
    double y = exp(-1.5*new_x_array[i]);
    double fHe = Yp/(4* (1-Yp));
    double yreion = pow((1+zreion), 1.5);
    double dyreion = 1.5 * sqrt(1+zreion) * dzreion;
    double intermediate = ((1+fHe)/2) * (1 + tanh((yreion-y)/dyreion));
    double Xe_corr = Xe_arr_peebles[i] + intermediate;

    // Append both the values with and without reionization to respective arrays
    Xe_arr.push_back(log(Xe_arr_peebles[i]));
    Xe_reion_arr.push_back(log(Xe_corr));
  };
  // Append the values with reionization correction to Xe_reion_arr
  
  

  log_Xe_of_x_spline.create(x_array, Xe_arr, "Function log_Xe(x)");
  log_Xe_reion_of_x_spline.create(x_array, Xe_reion_arr, "Function log_Xe(x) with reionization");
  log_Xe_of_x_Saha_spline.create(x_array, Xe_Saha_arr, "Function Xe(x) with Saha prediction");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB0 = cosmo->get_OmegaB(0);
  const double H0 = cosmo->get_H0();
  const double Tb = cosmo->get_TCMB(x);
  const double crit_density = (3*pow(H0,2))/(8*M_PI*G);

  // Compute baryon number density
  const double nb = (OmegaB0 * crit_density)/(m_H * pow(a,3.0));

  // Compute the R.H.S. of the Saha equation to apply the binary approximation if necessary
  const double R = 1/nb * pow((m_e*k_b*Tb)/(2*M_PI), 1.5) * pow(hbar, -3.0) * exp(-epsilon_0/(k_b*Tb));
  const double intermediate = (-R - sqrt(pow(R, 2.0) + 4.0*R));
  
  double Xe, ne;
  Xe = (-2*R)/intermediate;
  ne = nb * Xe;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double fine_structure = 1.0/137.0359992;

  const double OmegaB0 = cosmo->get_OmegaB(0);
  const double H0 = cosmo->get_H0();
  double Tb = cosmo->get_TCMB(x);
  double H = cosmo->H_of_x(x);
  const double crit_density = (3*pow(H0,2))/(8*M_PI*G);

  // Get number densities of hidrogen and 1s
  double nH = (1.0-Yp) * (OmegaB0 * crit_density)/(m_H * pow(a,3.0));
  double n1s = (1.0-X_e) * nH;

  // Compute Cr step by step (beta_2 is immediately defined without resorting to beta to avoid numerical errors)
  
  double phi_2 = 0.448 * log(epsilon_0/(k_b*Tb));
  double alpha_2 = (64.0*M_PI/sqrt(27.0*M_PI)) * ((pow(fine_structure, 2.0)*pow(hbar,2.0))/(c*pow(m_e,2.0))) * sqrt(epsilon_0/(k_b*Tb)) * phi_2;
  double beta = alpha_2 * pow((m_e*k_b*Tb)/(2.0*M_PI), 1.5) * pow(hbar, -3) * exp(-epsilon_0/(k_b*Tb));
  double beta_2 = alpha_2 * pow((m_e*k_b*Tb)/(2.0*M_PI), 1.5) * pow(hbar, -3) * exp(-epsilon_0/(4.0*k_b*Tb));
  double lambda_alpha = H * pow(3.0*epsilon_0, 3)/(pow(8.0*M_PI, 2.0)*n1s*pow(hbar*c, 3.0));
  double Cr = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta_2);
  
  // Fully define the R.H.S. of the equation
  dXedx[0] = (Cr/H) * (beta * (1.0-X_e) - nH * alpha_2 * pow(X_e,2.0));
  
  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){

  // Set up x-arrays to integrate over. Begin at x_end = 0 and integrate backwards
  const int npts = 1000;
  std::vector<Vector> x_arrays = {Utils::linspace(x_end, -log(1+zreion-dzreion), npts), Utils::linspace(-log(1+zreion-dzreion), -log(1+zreion+dzreion), npts_rec_arrays- 2 * npts),  Utils::linspace(-log(1+zreion+dzreion), x_start, npts)};


  // The ODE system dtau_noreion/dx, dtau/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    // Get free electron fractions and number densities in the 3 separate cases
    double Xe = Xe_of_x(x, false);
    double Xe_reion = Xe_of_x(x, true);
    double Xe_Saha = Xe_of_x_Saha(x);
    double ne = ne_of_x(x);
    double ne_reion = ne * (Xe_reion/Xe);
    double ne_Saha = ne * (Xe_Saha/Xe);
    double H = cosmo->H_of_x(x);

    
    //1st component: tau with no reionization
    //2nd component: tau with reionization
    //3rd component: tau only from Saha equation
    dtaudx[0] = -(Constants.c * ne * Constants.sigma_T)/H;
    dtaudx[1] = -(Constants.c * ne_reion * Constants.sigma_T)/H;
    dtaudx[2] = -(Constants.c * ne_Saha * Constants.sigma_T)/H;
    return GSL_SUCCESS;
  };

  //set up initial conditions
  Vector Tau_i = {0.0, 0.0, 0.0};

  //make vectors of vectors to solve the equation in each of the 3 x_arrays, then concatenate them
  std::vector<Vector> tau_arrays, tau_reion_arrays, tau_Saha_arrays;

  for (int i = 0; i<3; i++){
    ODESolver opticaldepth_ode;
    opticaldepth_ode.solve(dtaudx, x_arrays[i] , Tau_i, gsl_odeiv2_step_rk4);
    tau_arrays.push_back(opticaldepth_ode.get_data_by_component(0));
    tau_reion_arrays.push_back(opticaldepth_ode.get_data_by_component(1));
    tau_Saha_arrays.push_back(opticaldepth_ode.get_data_by_component(2));

    //redefine initial conditions for solving the equation in the next x_array
    Tau_i = {tau_arrays[i].back(), tau_reion_arrays[i].back(), tau_Saha_arrays[i].back()};
    // append the second and third x and tau arrays to the first (and begin the iteration from the 2nd element of each array)
    // to avoid duplicate x values
    if (i!=0){
    x_arrays[0].insert(x_arrays[0].end(), x_arrays[i].begin()+1, x_arrays[i].end()); 
    tau_arrays[0].insert(tau_arrays[0].end(), tau_arrays[i].begin()+1, tau_arrays[i].end());
    tau_reion_arrays[0].insert(tau_reion_arrays[0].end(), tau_reion_arrays[i].begin()+1, tau_reion_arrays[i].end());
    tau_Saha_arrays[0].insert(tau_Saha_arrays[0].end(), tau_Saha_arrays[i].begin()+1, tau_Saha_arrays[i].end());
    };
  };

  // reversing the arrays is necessary because otherwise splining does not work
  reverse(x_arrays[0].begin(), x_arrays[0].end());
  reverse(tau_arrays[0].begin(), tau_arrays[0].end());
  reverse(tau_reion_arrays[0].begin(), tau_reion_arrays[0].end());
  reverse(tau_Saha_arrays[0].begin(), tau_Saha_arrays[0].end());



  tau_of_x_spline.create(x_arrays[0], tau_arrays[0], "Function Tau(x)");
  tau_reion_of_x_spline.create(x_arrays[0], tau_reion_arrays[0], "Function Tau_reion(x)");
  tau_of_x_Saha_spline.create(x_arrays[0], tau_Saha_arrays[0], "Function Tau_Saha(x)");

  // Make the visibility function arrays and spline it (only for tau and tau_reion)
  Vector g_tilde_array, g_tilde_reion_array;

  auto visibility_function = [&] (const double x) {
    double tau = tau_of_x_spline(x);
    double tau_prime = tau_of_x_spline.deriv_x(x);

    double tau_reion = tau_reion_of_x_spline(x);
    double tau_reion_prime = tau_reion_of_x_spline.deriv_x(x);

    g_tilde_array.push_back(-tau_prime * exp(-tau));
    g_tilde_reion_array.push_back(-tau_reion_prime * exp(-tau_reion));
  };
  
  
  // Calculates the visibility function for all pts. in x_array and pushes it back to g_tilde_array
  for_each(x_arrays[0].begin(), x_arrays[0].end(), visibility_function);

  // Splines the result of the previous line
  g_tilde_of_x_spline.create(x_arrays[0], g_tilde_array, "Function g_tilde(x)");
  g_tilde_reion_of_x_spline.create(x_arrays[0], g_tilde_reion_array, "Function g_tilde_reion(x)");
}

//====================================================
// Solve for the sound horizon
//====================================================
void RecombinationHistory::solve_sound_horizon(){

  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);

  // Get the cosmological parameters needed for the calculation
  const double OmegaB0 = cosmo -> get_OmegaB(0);
  const double OmegaR0 = cosmo -> get_OmegaR(0);

  

  // The ODE for ds/dx
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
    // Get the scale factor and the conformal hubble factor, calculate the sound speed
    double a = exp(x);
    double Hp = cosmo -> Hp_of_x(x);
    double R = (4.0*OmegaR0)/(3.0*OmegaB0*a);
    double cs = Constants.c * sqrt(R/(3.0 * (1+R)));

    dsdx[0] = cs/Hp;
    return GSL_SUCCESS;
  };

  // Define the initial conditions
  const double ai = exp(x_start);
  double Hpi = cosmo -> Hp_of_x(x_start);
  double Ri = (4.0*OmegaR0)/(3.0*OmegaB0*ai);
  double csi = Constants.c * sqrt(Ri/(3.0 * (1+Ri)));
  Vector s_i = {csi/Hpi};

  // Solve the ODE
  ODESolver soundhorizon_ode;
  soundhorizon_ode.solve(dsdx, x_array , s_i, gsl_odeiv2_step_rk4);

  // Make sound horizon array and spline
  auto s_array = soundhorizon_ode.get_data_by_component(0);
  sound_horizon_spline.create(x_array, s_array, "Function s(x)");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x, bool reion) const{
  if(reion){
    return tau_reion_of_x_spline(x);
  } else{
    return tau_of_x_spline(x);
  }
}

double RecombinationHistory::dtaudx_of_x(double x, bool reion) const{
  if(reion){
    return tau_reion_of_x_spline.deriv_x(x);
  } else{
    return tau_of_x_spline.deriv_x(x);
  }
}

double RecombinationHistory::ddtauddx_of_x(double x, bool reion) const{
  if(reion){
    return tau_reion_of_x_spline.deriv_xx(x);
  } else{
    return tau_of_x_spline.deriv_xx(x);
  }
}

double RecombinationHistory::g_tilde_of_x(double x, bool reion) const{
  if(reion){
    return g_tilde_reion_of_x_spline(x);
  } else{
    return g_tilde_of_x_spline(x);
  }
  
}

double RecombinationHistory::dgdx_tilde_of_x(double x, bool reion) const{
  if(reion){
    return g_tilde_reion_of_x_spline.deriv_x(x);
  } else{
    return g_tilde_of_x_spline.deriv_x(x);
  }
}

double RecombinationHistory::ddgddx_tilde_of_x(double x, bool reion) const{
  if(reion){
    return g_tilde_reion_of_x_spline.deriv_xx(x);
  } else{
    return g_tilde_of_x_spline.deriv_xx(x);
  }
}

double RecombinationHistory::Xe_of_x(double x, bool reion) const{
  if(reion){
    return exp(log_Xe_reion_of_x_spline(x));
  }
  else{
    return exp(log_Xe_of_x_spline(x));
  }
}

double RecombinationHistory::Xe_of_x_Saha(double x) const{
  return exp(log_Xe_of_x_Saha_spline(x));
}

double RecombinationHistory::sound_horizon(double x) const{
  return sound_horizon_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  // Necessary cosmological parameters and scale factor
  const double a = exp(x);
  const double OmegaB0 = cosmo->get_OmegaB(0);
  const double H0 = cosmo->get_H0();
  const double crit_density = 3*pow(H0,2)/(8*M_PI*Constants.G);
  
  // Fetching Xe and determining baryon number density
  const double Xe = Xe_of_x(x, false);
  const double nb = (OmegaB0 * crit_density)/(Constants.m_H * pow(a,3));

  return Xe * nb;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

double RecombinationHistory::get_zreion() const{
  return zreion;
}

double RecombinationHistory::get_dzreion() const{
  return dzreion;
}

double RecombinationHistory::last_scattering_surface() const{
  return Utils::binary_search_for_value(tau_of_x_spline, 1.0, std::make_pair(x_start, x_end), 1e-8);
}

double RecombinationHistory::last_scattering_surface_Saha() const{
  return Utils::binary_search_for_value(tau_of_x_Saha_spline, 1.0, std::make_pair(x_start, x_end), 1e-8);
}

double RecombinationHistory::x_at_recombination() const{
  return Utils::binary_search_for_value(log_Xe_of_x_spline, log(0.1), std::make_pair(x_start, x_end), 1e-8);
}

double RecombinationHistory::x_at_recombination_Saha() const{
  return Utils::binary_search_for_value(log_Xe_of_x_Saha_spline, log(0.1), std::make_pair(x_start, x_end), 1e-8);
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << "Zreion:      " << zreion << "\n";
  std::cout << "dZreion:     " << dzreion << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = -15;
  const double x_max   = x_end;

  //Fetch L.S.S. and recombination both with full solution and only with Saha approximation
  double last_scattering = last_scattering_surface();
  double recombination = x_at_recombination();
  double last_scattering_Saha = last_scattering_surface_Saha();
  double recombination_Saha = x_at_recombination_Saha();

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  fp << "x   " << "Xe   " <<"ne     "<<"Tau     " <<"dTaudx   " <<"ddTauddx   "<<"g_tilde    "  <<"dgdx      "<<"ddgddx     "<<"Xe (Saha)  "<<"sound horizon"<<"(Reion) Xe " <<"Tau      " <<"dTaudx   " <<"ddTauddx   "<<"g_tilde    "  <<"dgdx      "<<"ddgddx     "<<std::endl;
  fp << std::endl;
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x, false)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x, false)          << " ";
    fp << dtaudx_of_x(x, false)       << " ";
    fp << ddtauddx_of_x(x, false)     << " ";
    fp << g_tilde_of_x(x, false)      << " ";
    fp << dgdx_tilde_of_x(x, false)   << " ";
    fp << ddgddx_tilde_of_x(x, false) << " ";
    fp << Xe_of_x_Saha(x) << " ";
    fp << sound_horizon(x) << " ";
    fp << Xe_of_x(x, true) << " ";
    fp << tau_of_x(x, true)          << " ";
    fp << dtaudx_of_x(x, true)       << " ";
    fp << ddtauddx_of_x(x, true)     << " ";
    fp << g_tilde_of_x(x, true)      << " ";
    fp << dgdx_tilde_of_x(x, true)   << " ";
    fp << ddgddx_tilde_of_x(x, true) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
  fp << "x          " << "t            " << "z              " << '\n';
  fp << "Last scattering surface:" << last_scattering << " "<< cosmo -> t_of_x(last_scattering) << " " << 1/exp(last_scattering)-1 << "\n";
  fp << "Recombination:" << recombination << " "<< cosmo -> t_of_x(recombination) << " " << 1/exp(recombination)-1 << " "<< cosmo->get_TCMB(recombination) <<"\n";
  fp << "LSS (Saha):" << last_scattering_Saha << " "<< cosmo -> t_of_x(last_scattering_Saha) << " " << 1/exp(last_scattering_Saha)-1 << "\n";
  fp << "Recombination (Saha):" << recombination_Saha << " "<< cosmo -> t_of_x(recombination_Saha) << " " << 1/exp(recombination_Saha)-1 << " "<< cosmo->get_TCMB(recombination_Saha) <<"\n";
  fp << "Sound horizon:" << sound_horizon(last_scattering);
}

