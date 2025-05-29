#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67; //0.67
  double OmegaB      = 0.05; //0.05
  double OmegaCDM    = 0.267; //0.267
  double OmegaK      = 0.0; //0.0
  double Neff        = 3.046; //3.046
  double TCMB        = 2.7255; //2.7255

  // Recombination parameters
  double Yp          = 0;
  double zreion = 8.0;
  double dzreion = 0.5;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965; //0.965
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  
  //mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp, zreion, dzreion);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue1 = 0.1 / Constants.Mpc, kvalue2 = 0.01 / Constants.Mpc, kvalue3 = 0.001 / Constants.Mpc;
  pert.output(kvalue1, "perturbations_k0.1.txt");
  pert.output(kvalue2, "perturbations_k0.01.txt");
  pert.output(kvalue3, "perturbations_k0.001.txt");
  
  
  //=========================================================================
  // Module IV
  //=========================================================================

  // Last variable corresponds to whether one wants to split the components of ST or not 
  // This takes quite a long time and is not strictly necessary, so set to zero if you do not want it
  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc, 0);
  power.solve();
  power.output("cells.txt", "pofk.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
