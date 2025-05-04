#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){
  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");
  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array(n_k);
  for (int ik = 0; ik < n_k; ik++) {
    k_array[ik] = pow(10, log10(k_min) + ik * (log10(k_max) - log10(k_min))/(n_k-1));
  }
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  const bool polarization = Constants.polarization;
  const bool neutrinos = Constants.neutrinos;

  // Constants and cosmological parameters
  const double c = Constants.c;
  const double OmegaR0 = cosmo -> get_OmegaR(0);
  const double OmegaB0 = cosmo -> get_OmegaB(0);
  const double OmegaCDM0 = cosmo -> get_OmegaCDM(0);
  const double OmegaNu0 = cosmo -> get_OmegaNu(0);
  const double H0 = cosmo -> get_H0();

  // Needed to implement reionization in case x > xreion
  const double zreion = rec -> get_zreion();
  const double xreion = log(1.0/(1.0+zreion));

  // Make a vector of vectors to contain the values of all evolved perturbations (n_ell_tot_full)
  // but also Pi and Psi (+ 2)
  std::vector<Vector> perturbation_arrays(Constants.n_ell_tot_full + 2);
  
  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Create the vectors that store the perturbations both in the tight coupling and full regimes
    // They are both set with n_ell_tot_full components to facilitate appending them later
    std::vector<Vector> y_tc(Constants.n_ell_tot_full), y_full(Constants.n_ell_tot_full);

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];
    // Find value to integrate to, split integration array into two intervals
    double x_end_tight = get_tight_coupling_time(k);

    Vector x_tc, x_full;

    
    bool tight_coupling_regime = 1;

    // For each value of x, check if we are still within the tight coupling regime
    for(int i=0; i<n_x; i++){
      if (x_array[i]<x_end_tight){
        x_tc.push_back(x_array[i]);

      } else {
          if(tight_coupling_regime){

            // When we leave tight coupling, append the x value at which it ends to both x arrays
            // (for the sake of consistency) and set the boolean to 0
            x_tc.push_back(x_end_tight);
            x_full.push_back(x_end_tight);
            tight_coupling_regime = 0;
          }
        x_full.push_back(x_array[i]);
      };
    }
  

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    ODESolver tight_coupling_sol;
    tight_coupling_sol.solve(dydx_tight_coupling, x_tc , y_tight_coupling_ini, gsl_odeiv2_step_rk4);

    // Fetch the solutions in the tight coupling regime component by component
    for (int i = 0; i<Constants.n_ell_tot_full; i++){

      if(i<Constants.ind_start_nu_tc){
        // If we haven't yet reached the neutrino multipoles, the component of the tight coupling solution
        // lines up with the component of the full system 
        y_tc[i] = tight_coupling_sol.get_data_by_component(i);

        // Between Constants.ind_start_nu_tc and Constants.ind_start_nu, append nothing (since these variables)
        // are only dynamically evolved in the full system

      } else if (i>=Constants.ind_start_nu && neutrinos){
        // If we have reached the neutrino multipoles, the component of the tight coupling solution
        // will differ from the one in the full system by Constants.n_ell_tot_full - Constants.n_ell_tot_tc
        y_tc[i] = tight_coupling_sol.get_data_by_component(i - (Constants.n_ell_tot_full - Constants.n_ell_tot_tc));
      }
      
    }
    

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    Vector y_tight_coupling;

    
    for(int i = 0; i<Constants.n_ell_tot_full; i++){
      if(i<Constants.ind_start_nu_tc || i>=Constants.ind_start_nu){ 
        // Append the last value of all y_tc arrays that contain elements
        // (smaller than ind_start_nu_tc and greater or equal than Constants.ind_start_nu)   
        y_tight_coupling.push_back(y_tc[i].back());
      }
    }

    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    ODESolver full_sol;
    full_sol.solve(dydx_full, x_full, y_full_ini);

    // Fetch the solutions in the full regime component by component
    for (int i = 0; i<Constants.n_ell_tot_full; i++){
      y_full[i] = full_sol.get_data_by_component(i);
    }

    //Append the arrays from the 2 regimes

    for (int i = 0; i<Constants.n_ell_tot_full; i++){

      // For the quantities that are dynamically evolved in both regimes, directly append
      // the y_full values to y_tc
      if(i<Constants.ind_start_nu_tc || i>=Constants.ind_start_nu){
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }

        
      // Photon quadrupole
      } else if (i==Constants.ind_start_theta + 2){
        
        // Append the values of Theta2 during tight coupling (it is not dynamically evolved in this regime)
        for(int j = 0; j<x_tc.size(); j++){

          // Background and recombination variables
          double Hp = cosmo -> Hp_of_x(x_tc[j]);
          double tauprime;

          // Calculate tauprime with a different function of x depending on whether we are past reionization or not
          if(x_tc[j] < xreion){
            tauprime = rec -> dtaudx_of_x(x_tc[j]);
          } else{
            tauprime = rec -> dtaudx_reion_of_x(x_tc[j]);
          }

          // For the tight coupling regime, simply append the value of Theta2 as given by the initial conditions to
          // y_tc[i] (this depends on whether or not polarization is true)
          if(polarization){
            y_tc[i].push_back(-(8.0*c*k)/(15.0*Hp*tauprime) * y_tc[Constants.ind_start_theta_tc + 1][j]);
          } else{
            y_tc[i].push_back(-(20.0*c*k)/(45.0*Hp*tauprime) * y_tc[Constants.ind_start_theta_tc + 1][j]);
          }
          
          
        }

        // Afterwards, directly append the values present in y_full, given by the ODESolver
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }

      // Higher order photon multipoles
      } else if ((i>Constants.ind_start_theta + 2) && (i<Constants.ind_start_thetap)){

        // Append the values of ThetaL during tight coupling (it is not dynamically evolved in this regime)
        for(int j = 0; j<x_tc.size(); j++){ 
          int l = i - Constants.ind_start_theta;
          double Hp = cosmo -> Hp_of_x(x_tc[j]);
          double tauprime;
          if(x_tc[j] < xreion){
            tauprime = rec -> dtaudx_of_x(x_tc[j]);
          } else{
            tauprime = rec -> dtaudx_reion_of_x(x_tc[j]);
          }
          y_tc[i].push_back(-(l)/(2.0*l + 1.0) * (c*k)/(Hp*tauprime) * y_tc[i-1][j]);
        }

        // Afterwards, directly append the values present in y_full, given by the ODESolver
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }
      
      // Polarization monopole
      } else if (i == Constants.ind_start_thetap){

        // Append the values of Thetap0 during tight coupling (it is not dynamically evolved in this regime)
        for(int j = 0; j < x_tc.size(); j++){
          y_tc[i].push_back(1.25 * y_tc[Constants.ind_start_theta + 2][j]);
        }

        // Afterwards, directly append the values present in y_full, given by the ODESolver
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }
      
      // Polarization dipole
      } else if (i == Constants.ind_start_thetap + 1){

        // Append the values of Thetap1 during tight coupling (it is not dynamically evolved in this regime)
        for(int j = 0; j < x_tc.size(); j++){
          double Hp = cosmo -> Hp_of_x(x_tc[j]);
          double tauprime;
          if(x_tc[j] < xreion){
            tauprime = rec -> dtaudx_of_x(x_tc[j]);
          } else{
            tauprime = rec -> dtaudx_reion_of_x(x_tc[j]);
          }
          y_tc[i].push_back((-c*k/(4*Hp*tauprime))*y_tc[Constants.ind_start_theta + 2][j]);
        }
        
        // Afterwards, directly append the values present in y_full, given by the ODESolver
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }
      
      // Polarization quadrupole
      } else if (i == Constants.ind_start_thetap + 2){

        // Append the values of Thetap2 during tight coupling (it is not dynamically evolved in this regime)
        for(int j = 0; j < x_tc.size(); j++){
          y_tc[i].push_back(0.25 * y_tc[Constants.ind_start_theta + 2][j]);
        }

        // Afterwards, directly append the values present in y_full, given by the ODESolver
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }
      
      // Higher order polarization multipoles
      } else if (i > Constants.ind_start_thetap + 2 && i < Constants.ind_start_nu) {

        int l = i - Constants.ind_start_thetap;

        // Append the values of ThetapL during tight coupling (it is not dynamically evolved in this regime)
        for(int j = 0; j < x_tc.size(); j++){
          double Hp = cosmo -> Hp_of_x(x_tc[j]);
          double tauprime;
          if(x_tc[j] < xreion){
            tauprime = rec -> dtaudx_of_x(x_tc[j]);
          } else{
            tauprime = rec -> dtaudx_reion_of_x(x_tc[j]);
          }
          y_tc[i].push_back(-l/(2.0*l + 1.0)* (c*k)/(Hp*tauprime) * y_tc[i-1][j]);
        }

        // Afterwards, directly append the values present in y_full, given by the ODESolver
        for(int j = 0; j<y_full[0].size(); j++){
          y_tc[i].push_back(y_full[i][j]);
        }

      }
    }
    

    // Append the values found in y_tc for this k value to the overall perturbation_arrays,
    // which encompasses all values of k

    for (int i=0; i<Constants.n_ell_tot_full+1; i++){

      for(int ix = 0; ix < n_x; ix++){

        if(i<Constants.n_ell_tot_full){
          // If i is still within the range of y_tc, simply append the values directly to 
          // perturbation_arrays[i]
          perturbation_arrays[i].push_back(y_tc[i][ix]);

        } else{
          // If past the range of y_tc, compute two extra values to append to perturbation_arrays (Pi and Psi)
          double a = exp(x_array[ix]);

          // If polarization is true, Pi also includes Thetap0 and Thetap2
          if(polarization){
            perturbation_arrays[i].push_back(y_tc[Constants.ind_start_theta + 2][ix] + y_tc[Constants.ind_start_thetap][ix] + y_tc[Constants.ind_start_thetap+2][ix]);
          } else {
            perturbation_arrays[i].push_back(y_tc[Constants.ind_start_theta + 2][ix]);
          }

          // If neutrinos is true, Psi also includes OmegaNu0 * N2
          double Psi = (neutrinos)? -y_tc[4][ix]-12.0*pow(H0, 2)/pow(c*k*a, 2)*(OmegaR0*y_tc[Constants.ind_start_theta+2][ix] + OmegaNu0*y_tc[Constants.ind_start_nu+2][ix]) : -y_tc[4][ix]-12.0*pow(H0, 2)/pow(c*k*a, 2)*OmegaR0*y_tc[Constants.ind_start_theta+2][ix];
          perturbation_arrays[i+1].push_back(Psi);
        }

      }
    }

  }
  Utils::EndTiming("integrateperturbation");

  // Create splines for the density, velocity and metric perturbations
  delta_cdm_spline.create(x_array, k_array, perturbation_arrays[0], "delta_cdm(x,k)");
  delta_b_spline.create(x_array, k_array, perturbation_arrays[1], "delta_b(x,k)");
  v_cdm_spline.create(x_array, k_array, perturbation_arrays[2], "v_cdm(x,k)");
  v_b_spline.create(x_array, k_array, perturbation_arrays[3], "v_b(x,k)");
  Phi_spline.create(x_array, k_array, perturbation_arrays[4], "Phi(x,k)");

  // The multipole splines have to be allocated before being used
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for(int i=0; i<Constants.n_ell_theta; i++){
    Theta_spline[i].create(x_array, k_array, perturbation_arrays[Constants.ind_start_theta+i],"Theta(x,k)");
  };

  // Only create the polarization and neutrino splines if their respective booleans are true
  if(polarization){
    Theta_p_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
    for(int i=0; i<Constants.n_ell_thetap; i++){
      Theta_p_spline[i].create(x_array, k_array, perturbation_arrays[Constants.ind_start_thetap+i],"Thetap(x,k)");
    };
  }

  if(neutrinos){
    Nu_spline = std::vector<Spline2D>(Constants.n_ell_neutrinos);
    for(int i=0; i<Constants.n_ell_neutrinos; i++){
      Nu_spline[i].create(x_array, k_array, perturbation_arrays[Constants.ind_start_nu+i],"Nu(x,k)");
    };
  }

  // Create two additional splines for the quantities Pi, Psi, which were not dynamically evolved
  Pi_spline.create(x_array, k_array, perturbation_arrays[Constants.n_ell_tot_full], "Pi(x,k)");
  Psi_spline.create(x_array, k_array, perturbation_arrays.back(), "Psi(x,k)");
  
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  
  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc.at(Constants.ind_deltacdm_tc);
  double &delta_b      =  y_tc.at(Constants.ind_deltab_tc);
  double &v_cdm        =  y_tc.at(Constants.ind_vcdm_tc);
  double &v_b          =  y_tc.at(Constants.ind_vb_tc);
  double &Phi          =  y_tc.at(Constants.ind_Phi_tc);
  double *Theta        = &y_tc.at(Constants.ind_start_theta_tc);
  
  // Constants and cosmological parameters
  const double c = Constants.c;
  const double a = exp(x);
  const double Hp = cosmo -> Hp_of_x(x);
  const double H0 = cosmo -> get_H0();
  const double OmegaR0 = cosmo -> get_OmegaR(0);
  const double OmegaNu0 = cosmo -> get_OmegaNu(0);

  // Define fv differently depending on whether neutrinos are included or not
  const double fv = (neutrinos) ? OmegaNu0/(OmegaR0 + OmegaNu0) : 0.0;


  double Psi = -(1.0)/(1.5 + 0.4*fv);
  Phi = -(1.0 + 0.4*fv)*Psi;
  delta_cdm = -1.5 * Psi;
  delta_b = delta_cdm;
  v_cdm = -((c*k)/(2*Hp)) * Psi;
  v_b = v_cdm;

  *Theta = -0.5 * Psi;
  *(Theta + 1) = (c*k)/(6*Hp) * Psi;

  // If neutrinos = true, set the initial conditions for neutrinos
  if(neutrinos){
    double *Nu = &y_tc.at(Constants.ind_start_nu_tc);

    *Nu = -0.5*Psi;
    *(Nu + 1) = (c*k)/(6*Hp) * Psi;
    *(Nu + 2) = -pow(c*k*a, 2)*(Psi + Phi)/(12 * pow(H0, 2)*OmegaNu0);

    for (int l = 3; l<n_ell_neutrinos_tc; l++){
      *(Nu + l) = (c*k)/((2.0*l + 1.0) * Hp) * *(Nu + l - 1);
    }
  } 
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);

  // A double to point to with the Nu, dNudx, Thetap, dThetapdx components if neutrinos = false and/or polarization = false
  double null = 0.0;

  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc.at(Constants.ind_deltacdm_tc);
  const double &delta_b_tc      =  y_tc.at(Constants.ind_deltab_tc);
  const double &v_cdm_tc        =  y_tc.at(Constants.ind_vcdm_tc);
  const double &v_b_tc          =  y_tc.at(Constants.ind_vb_tc);
  const double &Phi_tc          =  y_tc.at(Constants.ind_Phi_tc);
  const double *Theta_tc        = &y_tc.at(Constants.ind_start_theta_tc);
  

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Thetap = (polarization)? &y[Constants.ind_start_thetap] : &null;
  

  
  
  const double c = Constants.c;
  
  // Needed to implement reionization in case x > xreion
  const double zreion = rec -> get_zreion();
  const double xreion = log(1.0/(1.0+zreion));
  double tauprime;
  if(x < xreion){
    tauprime = rec -> dtaudx_of_x(x);
  } else{
    tauprime = rec -> dtaudx_reion_of_x(x);
  }

  // Background variables
  const double Hp = cosmo -> Hp_of_x(x);
  const double eta = cosmo -> eta_of_x(x);

  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;
  *Theta = *Theta_tc;
  *(Theta + 1) = *(Theta_tc + 1);

  // If polarization is true, set the initial conditions for the polarization multipoles
  // The initial value of the photon quadrupole also depends on whether polarization is true or false
  if(polarization){
    *(Theta + 2) = -(8.0*c*k)/(15.0*Hp*tauprime) * *(Theta_tc + 1);

    *Thetap = 1.25 * *(Theta + 2);
    *(Thetap + 1) = -(c*k)/(4.0*Hp*tauprime) * *(Theta + 2);
    *(Thetap + 2) = 0.25* *(Theta + 2);
    for (int l = 3; l<n_ell_thetap; l++){
      *(Thetap + l) = (-l/(2.0*l + 1.0)) * (c*k)/(Hp*tauprime) * *(Thetap + l-1);
    }
  } else{
    *(Theta + 2) = -(20.0*c*k)/(45.0*Hp*tauprime) * *(Theta_tc + 1);
  }

  for (int l = 3; l<n_ell_theta; l++){
    *(Theta + l) = (-l/(2.0*l+1.0)) * (c*k)/(Hp*tauprime) * *(Theta + l-1);
  }

  // If neutrinos is true, set the initial conditions for the neutrino multipoles
  if(neutrinos){
    const double *Nu_tc = &y_tc.at(Constants.ind_start_nu_tc);
    double *Nu = &y[Constants.ind_start_nu];
    for (int l = 0; l<n_ell_neutrinos; l++){
      *(Nu + l) = *(Nu_tc+l);
    }
  }



  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{

  // Variable declarations and constants
  double x_tight_coupling_end;
  double tauprime; 
  const double c = Constants.c;

  // Needed to implement reionization in case x > xreion
  const double zreion = rec -> get_zreion();
  const double xreion = log(1.0/(1.0+zreion));
  
  Spline tight_coupling_spline; 
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  Vector tight_coupling_array(n_x);
  for(size_t i = 0; i < x_array.size(); i++){
    // Compute the derivative of the optical depth depending on if we are past reionization or not
    if(x_array[i]< xreion){
      tauprime = rec -> dtaudx_of_x(x_array[i]);
    } else{
      tauprime = rec -> dtaudx_reion_of_x(x_array[i]);
    }
    
    double Hp = cosmo -> Hp_of_x(x_array[i]);

    // Check which of c*k/Hp or 1 is bigger
    if ((c*k/Hp) < 1){
      // Subtract the absolute value of tauprime from 10 times max(c*k/Hp, 1)
      tight_coupling_array[i] = fabs(tauprime) - 10.0;
    } else{
      tight_coupling_array[i] = fabs(tauprime) - 10.0 * (c*k/Hp);
    };
  }

  // Create a spline of |tauprime|-10*max(c*k/Hp, 1), search for 0 to find the end of tight coupling
  tight_coupling_spline.create(x_array, tight_coupling_array, "Spline to test tight coupling");

  double x_tight_coupling_test = Utils::binary_search_for_value(tight_coupling_spline, 0.0, std::make_pair(x_start, x_end), 1e-8);

  // If we have already left recombination, end tight coupling at recombination instead
  if (x_tight_coupling_test > -8.3){
    x_tight_coupling_end = -8.3;
  } else{
    x_tight_coupling_end = x_tight_coupling_test;
  }

  return x_tight_coupling_end;
}

//====================================================
// After integrating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");
  Vector k_array(n_k);
  for(int ik=0; ik < n_k; ik++){
    k_array[ik] = k_min + ik * (k_max-k_min)/(n_k-1);
  };

  Vector x_array=Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Necessary constants
  const double c = Constants.c;
  const double eta0 = cosmo -> eta_of_x(0);

  // Needed to implement reionization in case x > xreion
  const double zreion = rec -> get_zreion();
  const double xreion = log(1.0/(1.0+zreion));

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //Fetching background cosmology values
      const double Hp = cosmo->Hp_of_x(x);
      const double Hpp = cosmo -> dHpdx_of_x(x);
      const double Hppp = cosmo -> ddHpddx_of_x(x);
      const double eta = cosmo -> eta_of_x(x);
      
      //Fetching recombination values depending on if we are past reionization or not
      const double tau = (x < xreion)? rec->tau_of_x(x) : rec->tau_reion_of_x(x);
      const double g_tilde = (x < xreion)? rec->g_tilde_of_x(x) : rec->g_tilde_reion_of_x(x);
      const double g_tilde_prime = (x < xreion)? rec->dgdx_tilde_of_x(x) : rec->dgdx_tilde_reion_of_x(x);
      const double g_tilde_pp = (x < xreion)? rec->ddgddx_tilde_of_x(x) : rec->ddgddx_tilde_reion_of_x(x);

      //Getting all necessary values from the perturbation splines
      const double psi = get_Psi(x,k);
      const double psi_prime = Psi_spline.deriv_x(x,k);
      const double phi_prime = Phi_spline.deriv_x(x,k);
      const double pi = get_Pi(x,k);
      const double pi_prime = Pi_spline.deriv_x(x,k);
      const double pi_pp = Pi_spline.deriv_xx(x,k);
      const double v_b = get_v_b(x,k);
      const double v_b_prime = v_b_spline.deriv_x(x,k);
      const double theta_0 = get_Theta(x,k,0);

      // Simple application of the product rule to compute ddx(Hp g_tilde v_b) and ddx(Hp ddx(Hp g_tilde pi))
      const double dHp_g_v = Hpp * g_tilde * v_b + Hp * g_tilde_prime * v_b + Hp * g_tilde * v_b_prime;
      const double ddHp_g_pi = pow(Hpp, 2)*g_tilde*pi + Hppp*Hp*g_tilde*pi + 3* Hpp * Hp * (g_tilde_prime*pi + g_tilde*pi_prime) + pow(Hp, 2) * (2*g_tilde_prime*pi_prime + g_tilde_pp*pi + g_tilde*pi_pp);



      // Temperatur source
      ST_array[index] = g_tilde * (theta_0 + psi + (pi/4.0)) + exp(-tau)*(psi_prime-phi_prime)- dHp_g_v/(c*k)+ 3.0/(4.0*pow(c*k, 2))*ddHp_g_pi;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = (3 * g_tilde * pi)/(4 * pow(k, 2) * pow(eta0 - eta, 2));
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;
  const bool polarization          = Constants.polarization;

  // A double to point to with the Nu and dNudx components if neutrinos = false
  double null = 0.0;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu = (neutrinos)? &y[Constants.ind_start_nu_tc] : &null;

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx = (neutrinos)? &dydx[Constants.ind_start_nu_tc] : &null;

  // Constants and background variables
  const double c = Constants.c;
  const double OmegaR0 = cosmo -> get_OmegaR(0);
  const double OmegaB0 = cosmo -> get_OmegaB(0);
  const double OmegaCDM0 = cosmo -> get_OmegaCDM(0);
  const double OmegaNu0 = cosmo -> get_OmegaNu(0); 
  const double H0 = cosmo -> get_H0();
  double a = exp(x);
  double eta = cosmo->eta_of_x(x);
  const double R = (4.0*OmegaR0)/(3.0*OmegaB0*a);
  
  const double Hp = cosmo -> Hp_of_x(x);
  const double Hpp = cosmo -> dHpdx_of_x(x);
  

  // Define first two derivatives of tau depending on whether reionization has happened
  double tauprime;
  double taupprime;
  const double zreion = rec -> get_zreion();
  const double xreion = log(1.0/(1.0+zreion));

  if(x < xreion){
    tauprime = rec -> dtaudx_of_x(x);
    taupprime = rec -> ddtauddx_of_x(x);
  } else{
    tauprime = rec -> dtaudx_reion_of_x(x);
    taupprime = rec -> ddtauddx_reion_of_x(x);
  }


  ////////////////////////////////////////////////////////////// 
  // Additional definitions that are not dinamically evolved
  /////////////////////////////////////////////////////////////

  // Theta2 is different depending on whether polarization is true
  double Theta2 = (polarization)? -(8.0*c*k)/(15.0*Hp*tauprime) * *(Theta + 1) : -(20*c*k)/(45*Hp*tauprime) * *(Theta + 1);

  // If neutrinos is true, Psi and intermediate also include OmegaNu0 * N2
  double Psi = (neutrinos)? -Phi - (12*pow(H0, 2))/pow(c*k*a, 2) * (OmegaR0 * Theta2 + OmegaNu0 * *(Nu + 2)) : -Phi - (12*pow(H0, 2))/pow(c*k*a, 2) * (OmegaR0 * Theta2);
  double intermediate = (neutrinos)? OmegaCDM0*delta_cdm*pow(a, -1) + OmegaB0*delta_b*pow(a, -1) + 4*OmegaR0*(*Theta)*pow(a, -2) + 4*OmegaNu0*(*Nu)*pow(a, -2) : OmegaCDM0*delta_cdm*pow(a, -1) + OmegaB0*delta_b*pow(a, -1) + 4*OmegaR0*(*Theta)*pow(a, -2);
  

  dPhidx= Psi - pow(c*k, 2)/(3*pow(Hp, 2))*Phi + (pow(H0, 2)/(2*pow(Hp, 2)))*intermediate;
  *dThetadx = -(c*k/Hp)* *(Theta + 1) - dPhidx;

  // Define q to insert in the equation of dv_bdx
  double q = (-(3.0 * *(Theta+1) + v_b)*((1.0-R)*tauprime + (1.0+R)*taupprime) - (c*k/Hp)*Psi + (1.0 - Hpp/Hp)*(c*k/Hp)*(-*Theta + 2.0*Theta2) - c*k/Hp * *dThetadx)/((1.0+R)*tauprime + (Hpp/Hp) - 1.0);

  ddelta_cdmdx = (c*k/Hp) * v_cdm - 3.0*dPhidx;
  ddelta_bdx = (c*k/Hp) * v_b - 3.0*dPhidx;
  dv_cdmdx = -v_cdm - (c*k/Hp) * Psi;
  dv_bdx= (1.0/(1.0+R))*(-v_b - (c*k/Hp)*Psi + R * (q + (c*k/Hp)*(-*Theta + 2.0*Theta2) - (c*k/Hp)*Psi));

  *(dThetadx + 1)= (1.0/3.0) * (q - dv_bdx);

  // If neutrinos is true, d
  if(neutrinos){
    *dNudx = -(c*k/Hp)* *(Nu + 1) - dPhidx;
    *(dNudx + 1) = (c*k)/(3*Hp)* *Nu - (2*c*k)/(3*Hp)* *(Nu + 2) + (c*k)/(3*Hp) * Psi;
    for (int l = 2; l<n_ell_neutrinos_tc-1 ; l++){
      *(dNudx + l) = (l*c*k)/((2.0*l+1.0)*Hp)* *(Nu + l-1) - ((l+1.0)*c*k)/((2.0*l+1.0)*Hp)* *(Nu + l+1);
    }

    *(dNudx + n_ell_neutrinos_tc-1) = (c*k/Hp)* *(Nu + n_ell_neutrinos_tc-2) - c*(n_ell_neutrinos_tc)/(Hp*eta)* *(Nu + n_ell_neutrinos_tc-1);
  }
  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;


  // A double to point to with the Nu, dNudx, Thetap, dThetapdx components if neutrinos = false and/or polarization = false
  double null = 0.0;
  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Nu              =  (neutrinos)? &y[Constants.ind_start_nu] : &null;
  const double *Thetap          = (polarization)? &y[Constants.ind_start_thetap] : &null;

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dNudx           = (neutrinos)? &dydx[Constants.ind_start_nu] : &null;
  double *dThetapdx       = (polarization)? &dydx[Constants.ind_start_thetap]: &null;

  
  // Cosmological parameters and variables
  double a = exp(x);
  double Hp = cosmo->Hp_of_x(x);
  double eta = cosmo->eta_of_x(x);
  const double H0 = cosmo->get_H0();
  const double c = Constants.c;
  const double OmegaCDM0 = cosmo->get_OmegaCDM(0);
  const double OmegaB0 = cosmo->get_OmegaB(0);
  const double OmegaR0 = cosmo->get_OmegaR(0);
  const double OmegaNu0 = cosmo->get_OmegaNu(0);
  const double R = (4.0*OmegaR0)/(3.0*OmegaB0*a);

  // Recombination variables
  double tauprime;
  const double zreion = rec -> get_zreion();
  const double xreion = log(1.0/(1.0+zreion));
  // 
  if(x < xreion){
    tauprime = rec -> dtaudx_of_x(x);
  } else{
    tauprime = rec -> dtaudx_reion_of_x(x);
  }


  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi = (neutrinos)? -Phi - (12*pow(H0, 2))/pow(c*k*a, 2) * (OmegaR0 * *(Theta + 2) + OmegaNu0 * *(Nu + 2)) : -Phi - (12*pow(H0, 2))/pow(c*k*a, 2) * (OmegaR0 * *(Theta + 2));
  double Pi = (polarization)? *(Theta + 2) + *Thetap + *(Thetap + 2) : *(Theta + 2);
  double intermediate = (neutrinos)? OmegaCDM0*delta_cdm*pow(a, -1) + OmegaB0*delta_b*pow(a, -1) + 4*OmegaR0*(*Theta)*pow(a, -2) + 4*OmegaNu0*(*Nu)*pow(a, -2) : OmegaCDM0*delta_cdm*pow(a, -1) + OmegaB0*delta_b*pow(a, -1) + 4*OmegaR0*(*Theta)*pow(a, -2);
  
  dPhidx = Psi - pow(c*k, 2)/(3.0*pow(Hp, 2))*Phi + pow(H0, 2)/(2.0*pow(Hp, 2))*intermediate;
  ddelta_cdmdx = (c*k/Hp) * v_cdm - 3.0*dPhidx;
  ddelta_bdx = (c*k/Hp) * v_b - 3.0*dPhidx;
  dv_cdmdx = -v_cdm - (c*k/Hp) * Psi;
  dv_bdx= -v_b - (c*k/Hp) * Psi + tauprime*R*(3* *(Theta + 1) + v_b);

  *dThetadx = -(c*k/Hp)* *(Theta + 1) - dPhidx;
  *(dThetadx + 1) = (c*k)/(3.0*Hp)* *(Theta) - (2*c*k)/(3.0*Hp)* *(Theta + 2) + (c*k)/(3.0*Hp)*Psi + tauprime*(*(Theta + 1) + v_b/3.0);
  for (int l = 2; l<n_ell_theta-1 ; l++){
    if (l==2){
      *(dThetadx + l) = (l*c*k)/((2.0*l+1.0)*Hp)* *(Theta + l-1) - ((l+1.0)*c*k)/((2.0*l+1.0)*Hp)* *(Theta + l+1) + tauprime*(*(Theta + l) + 0.1*Pi);
    } else{
      *(dThetadx + l) = (l*c*k)/((2.0*l+1.0)*Hp)* *(Theta + l-1) - ((l+1.0)*c*k)/((2.0*l+1.0)*Hp)* *(Theta + l+1) + tauprime*(*(Theta + l));
    }
  };
  *(dThetadx + n_ell_theta-1) = (c*k/Hp)* *(Theta + n_ell_theta-2) - c*(n_ell_theta)/(Hp*eta)* *(Theta + n_ell_theta-1) + tauprime* *(Theta+n_ell_theta-1);


  if(neutrinos){
    *dNudx = -(c*k/Hp)* *(Nu + 1) - dPhidx;
    *(dNudx + 1) = (c*k)/(3*Hp)* *Nu - (2*c*k)/(3*Hp)* *(Nu + 2) + (c*k)/(3*Hp) * Psi;
    for (int l = 2; l<n_ell_neutrinos-1 ; l++){
      *(dNudx + l) = (l*c*k)/((2.0*l+1.0)*Hp)* *(Nu + l-1) - ((l+1.0)*c*k)/((2.0*l+1.0)*Hp)* *(Nu + l+1);
    }

    *(dNudx + n_ell_neutrinos-1) = (c*k/Hp)* *(Nu + n_ell_neutrinos-2) - c*(n_ell_neutrinos)/(Hp*eta)* *(Nu + n_ell_neutrinos-1);
  }

  if(polarization){
    *dThetapdx = -(c*k/Hp)* *(Thetap + 1) + tauprime * (*Thetap - 0.5*Pi);
    for (int l = 1; l<n_ell_thetap-1 ; l++){
      if (l==2){
        *(dThetapdx + l) = (l*c*k)/((2.0*l+1.0)*Hp)* *(Thetap + l-1) - ((l+1.0)*c*k)/((2.0*l+1.0)*Hp)* *(Thetap + l+1) + tauprime * (*(Thetap + l) - 0.1*Pi);
      } else{
        *(dThetapdx + l) = (l*c*k)/((2.0*l+1.0)*Hp)* *(Thetap + l-1) - ((l+1.0)*c*k)/((2.0*l+1.0)*Hp)* *(Thetap + l+1) + tauprime * *(Thetap + l); 
      }
    }
    *(dThetapdx + n_ell_thetap-1) = (c*k/Hp)* *(Thetap + n_ell_thetap-2) - c*(n_ell_thetap)/(Hp*eta)* *(Thetap + n_ell_thetap-1) + tauprime* *(Thetap + n_ell_thetap-1);
  }
  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;

  // Fetch times of matter radiation and matter dark energy equality, last scattering surface
  double matradequal = cosmo->matter_radiation_equality();
  double matlambdaequal = cosmo->matter_darkenergy_equality();
  double last_scattering = rec->last_scattering_surface();

  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Theta_p(x,k,0)   << " ";
    fp << get_Theta_p(x,k,1)   << " ";
    fp << get_Theta_p(x,k,2)   << " ";
    fp << get_Nu(x,k,0)   << " ";
    fp << get_Nu(x,k,1)   << " ";
    fp << get_Nu(x,k,2)   << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << get_delta_b(x,k) << " ";
    fp << get_delta_cdm(x,k) << " ";
    fp << get_v_b(x,k) << " ";
    fp << get_v_cdm(x,k) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
  fp << "x          " << "t            " << "z              " << '\n';
  fp << "Matter radiation equality:" << matradequal << " " << cosmo->t_of_x(matradequal) <<" "<< 1/exp(matradequal)-1 << "\n";
  fp << "Matter dark energy equality:" << matlambdaequal << " " << cosmo->t_of_x(matlambdaequal) <<" "<< 1/exp(matlambdaequal)-1 << "\n";
  fp << "Last scattering surface:" << last_scattering << " "<< cosmo -> t_of_x(last_scattering) << " " << 1/exp(last_scattering)-1 << "\n";
}

