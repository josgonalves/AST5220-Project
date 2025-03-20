#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  H0 = h * Constants.H0_over_h;
  OmegaR = (2 * pow(M_PI,3) * pow(Constants.k_b*TCMB,4) * 8 * Constants.G)/(90 * pow(Constants.hbar,3) * pow(Constants.c,5) * pow(H0,2));
  OmegaNu = Neff * (7.0/8.0) * pow(4.0/11.0, 4.0/3.0) * OmegaR;
  OmegaLambda = 1 - (OmegaB + OmegaCDM + OmegaK + OmegaR + OmegaNu);
  
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  //Defining the integration start and end points
  const double xmin = Constants.x_start;
  const double xmax = Constants.x_end;
  const double npts = 500;

  Vector x_array = Utils::linspace(xmin, xmax, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    double Hp = Hp_of_x(x);

    detadx[0] = Constants.c/Hp;

    return GSL_SUCCESS;
  };

  //initial conditions
  double eta0 = Constants.c/Hp_of_x(xmin);
  Vector eta_i = {eta0};

  //solving of the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_i, gsl_odeiv2_step_rk4);

  //get the solution and make the spline
  auto eta_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "Function eta(x)");
  

  //Now solving for the cosmic time, t

  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    double H = H_of_x(x);

    dtdx[0] = 1/H;

    return GSL_SUCCESS;
  };

  //initial conditions
  double t0 = 1/(2*H_of_x(xmin));
  
  Vector t_i = {t0};

  //Solving of the ODE
  ODESolver ode2;
  ode2.solve(dtdx, x_array, t_i, gsl_odeiv2_step_rk4);

  //get the solution and make the spline
  auto t_array = ode2.get_data_by_component(0);
  t_of_x_spline.create(x_array, t_array, "Function t(x)");

}

//====================================================
// Get methods
//====================================================


double BackgroundCosmology::H_of_x(double x) const{
  //define scale factor from x 
  double a = exp(x);
  //contributions from the terms with different equations of state
  double BandCDM = (OmegaB + OmegaCDM) * pow(a,-3);
  double RandNu = (OmegaR + OmegaNu) * pow(a,-4);
  double Curvature = OmegaK * pow(a,-2);
  double Lambda = OmegaLambda;

  return H0 * sqrt(BandCDM + RandNu + Curvature + Lambda);
}

double BackgroundCosmology::Hp_of_x(double x) const{
  //define hubble parameter and scale factor from x (this will also show up in all subsequent methods)
  double H = H_of_x(x);
  double a = exp(x);

  return H * a;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  //derivative calculated by the product rule, i.e. H * dadx + a * dHdx
  double H = H_of_x(x);
  double a = exp(x);

  //derivative of the Hubble parameter with respect to x
  double intermediate = (-3* (OmegaB + OmegaCDM) * pow(a,-3) - 4 * (OmegaR + OmegaNu) * pow(a,-4) - 2 * OmegaK * pow(a,-2) );
  double dHdx = (pow(H0,2) * intermediate)/(2 * H);

  //follows from the product rule, as dadx = a
  return a * (H + dHdx);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  double H = H_of_x(x);
  double a = exp(x);

  //using the derivatives of both H and Hp with respect to x simplifies the expression greatly
  double intermediate = (-3* (OmegaB + OmegaCDM) * pow(a,-3) - 4 * (OmegaR + OmegaNu) * pow(a,-4) - 2 * OmegaK * pow(a,-2) );
  double dHdx = (pow(H0,2) * intermediate)/(2 * H);
  double dHpdx = dHpdx_of_x(x);

  //second derivative of the hubble parameter w.r.t. x
  double ddHddx = (pow(H0,2)/2) * ( (9* (OmegaB + OmegaCDM) * pow(a,-3) + 16 * (OmegaR + OmegaNu) * pow(a,-4) + 4 * OmegaK * pow(a,-2) )/H - (dHdx * intermediate)/pow(H,2));
  return dHpdx + a*(dHdx + ddHddx);
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  double a = exp(x);
  double H = H_of_x(x);

  return OmegaB * pow(a,-3) * pow(H0/H, 2);
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  double a = exp(x);
  double H = H_of_x(x);

  return OmegaR * pow(a,-4) * pow(H0/H, 2);
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  double a = exp(x);
  double H = H_of_x(x);

  return OmegaNu * pow(a,-4) * pow(H0/H, 2);
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  double a = exp(x);
  double H = H_of_x(x);

  return OmegaCDM * pow(a,-3) * pow(H0/H, 2);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double H = H_of_x(x);

  return OmegaLambda * pow(H0/H, 2);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  double a = exp(x);
  double H = H_of_x(x);

  return OmegaK * pow(a,-2) * pow(H0/H, 2);
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  double cdistance = get_comoving_distance_of_x(x);
  double r;

  if (OmegaK == 0.0){
    r = cdistance;
  }
  else if (OmegaK > 0.0){
    double intermediate = cdistance*sqrt(abs(OmegaK))*H0/Constants.c;
    r = cdistance * (sinh(intermediate)/intermediate); 
  }
  else if (OmegaK < 0.0){
    double intermediate = cdistance*sqrt(abs(OmegaK))*H0/Constants.c;
    r = cdistance * (sin(intermediate)/intermediate); 
  }

  double a = exp(x);

  return r*a;
}

double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{

  double adistance = get_angular_distance_of_x(x);
  double a = exp(x);

  return adistance/pow(a,2);
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{

  double eta0 = eta_of_x(0);
  double eta = eta_of_x(x);

  return eta0 - eta;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::eta_prime_of_x(double x) const{
  return eta_of_x_spline.deriv_x(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

// These are simply analytical expressions for the value of x for which
// OmegaM = OmegaRad. and OmegaM = OmegaLambda
double BackgroundCosmology::matter_radiation_equality() const{
  return log((OmegaR+OmegaNu)/(OmegaB+OmegaCDM));
}

double BackgroundCosmology::matter_darkenergy_equality() const{
  return log((OmegaB+OmegaCDM)/OmegaLambda)/3;
}

// The acceleration starts when the second time derivative of the
// scale factor becomes >0. This is equivalent to looking for
// the value of x for which dHpdx becomes > 0

double BackgroundCosmology::start_of_acceleration() const{
  Spline dHp_of_x_spline; 
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, 500);
  Vector dHp_of_x_array(500);

  for(size_t i = 0; i < x_array.size(); i++) 
    dHp_of_x_array[i] = dHpdx_of_x(x_array[i]);

  dHp_of_x_spline.create(x_array, dHp_of_x_array, "Function dHpdx(x)");
  
  return Utils::binary_search_for_value(dHp_of_x_spline, 0, std::make_pair(Constants.x_start, Constants.x_end), 1e-8);
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << "Sum:         " << OmegaB + OmegaK + OmegaCDM + OmegaLambda + OmegaNu + OmegaR << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -15.0;
  const double x_max =  5.0;
  const int    n_pts =  1500;
  
  double matradequal = matter_radiation_equality();
  double matlambdaequal = matter_darkenergy_equality();
  double accelstart = start_of_acceleration();

  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  fp << "x   " << "eta   " <<"t     " <<"Hp   " << "dHpdx   " << "ddHpddx    " <<"dL      "<<"OmegaB   " << "OmegaCDM    " 
  << "OmegaLambda    " << "OmegaR    " << "OmegaNu    "<< "OmegaK   " << "eta prime    " << std::endl;
  fp << std::endl;
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << t_of_x(x)          << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)    << " ";
    fp << get_luminosity_distance_of_x(x)<<" ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << eta_prime_of_x(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
  fp << "Current age of the Universe:"<< eta_of_x(0) << " " << t_of_x(0) << " " <<'\n';
  fp << "x          " << "t            " << "z              " << '\n';
  fp << "Matter radiation equality:" << matradequal << " " << t_of_x(matradequal) <<" "<< 1/exp(matradequal)-1 << "\n";
  fp << "Matter dark energy equality:" << matlambdaequal << " " << t_of_x(matlambdaequal) <<" "<< 1/exp(matlambdaequal)-1 << "\n";
  fp << "Start of acceleration:" << accelstart << " " << t_of_x(accelstart) <<" "<< 1/exp(accelstart)-1 << "\n";
} 

