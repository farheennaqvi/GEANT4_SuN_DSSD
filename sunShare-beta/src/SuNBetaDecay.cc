// Alex Dombos
// August 2015
// Purpose: To calculate all quantites related to beta decay
//          needed for GEANT simulations with SuN.
//
//          Journal articles or Textbooks referenced in this program:
//          ASR:
//             Title: A simple relation for the Fermi function
//             Author: P. Venkataramaiah et. al
//             Journal: J. Phys. G: Nucl. Phys.
//          RDG:
//             Title: Radioactive Decays in Geant4
//             Author: S. Hauf et. al
//             Journal: arXiv
//          BDP
//             Title: Beta-decay properties of neutron-rich Zr and Mo isotopes
//             Author: P. Sarriguren and J. Pereira
//             Journal: Phys Rev C 81, 064314 (2010)
//
//          Websites that were helpful for writing one of the approximate Fermi functions
//          http://hurel.hanyang.ac.kr/Geant4/Doxygen/10.00/html/_g4_beta_fermi_function_8cc_source.html
//          http://geant4.cern.ch/support/source/geant4/source/processes/hadronic/models/radioactive_decay/src/G4BetaFermiFunction.cc
//          http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classes.html
//          http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classG4BetaDecayCorrections.html
//          http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classG4BetaFermiFunction.html
//          http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classG4BetaMinusDecayChannel.html
//          http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classG4BetaPlusDecayChannel.html
//
//          Always perform some test cases. Perform some
//          simulations with some excited states and plot the
//          exact precise analytical function on top of the
//          distribution from GEANT to make sure they are identical.
//
//          Helpful G4cout statements are commented out.
//          If uncommented, they will help you make sure
//          the code is performing correctly.
//
//          The following variables are defined in SuNBetaDecay.hh:
//          zDN
//          aDN
//          Qvalue
//          electronRestMassEnergy
//          zDN, aDN, and Qvalue will need to be changed depending on the beta decay

#include "SuNBetaDecay.hh"
#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh" 
#include "TMath.h"
#include <complex>

// Calculate the Gamma function for complex argument a+bi
// From http://en.wikipedia.org/wiki/Lanczos_approximation
std::complex<G4double> SuNBetaDecay::ComplexGammaFunction(const G4double a, const G4double b){
  // Coefficients used by the GNU Scientific Library
  G4double g = 7.0;
  G4double p[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
		   771.32342877765313, -176.61502916214059, 12.507343278686905,
		   -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
  
  std::complex<double> zcn(a,b); // (c = complex, n = number)
  // Reflection formula
  if (zcn.real() < 0.5) {
    std::complex<double> tmp = 1.0 - zcn;
    return TMath::Pi() / (sin(TMath::Pi()*zcn) * ComplexGammaFunction(tmp.real(), tmp.imag()));
  }
  else {
    zcn -= 1.0;
    std::complex<double> x = p[0];
    for (int i=1; i<g+2; i++){
      x += p[i]/(zcn+double(i));
    }
    std::complex<double> t = zcn + g + 0.5;
    std::complex<double> cG = sqrt(2*TMath::Pi()) * pow(t, zcn+0.5) * exp(-t) * x;
    return cG;
  }
}

G4double SuNBetaDecay::GetMaximumElectronKineticEnergy(){
  G4double sumOfGammaRayEnergies = 0.0;
  for (int i=1; i<=nEnergies; i++){ // Sum the energies of the gamma rays for a given cascade in the input file. Start at i=1 because we need to skip over the intensity which is in the 0th slot of the inputArray
    sumOfGammaRayEnergies = sumOfGammaRayEnergies + inputArray[cascade][i];
    //G4cout << "Gamma Ray Energy: " << inputArray[cascade][i] << "   " << "Cumulative Sum of Gamma Ray Energies: " << sumOfGammaRayEnergies  << G4endl;
  }
  G4double excitedStateEnergy = sumOfGammaRayEnergies;
  G4double maxKE = (Qvalue - excitedStateEnergy) / electronRestMassEnergy; // Calculate the maximum kinetic energy available to the electron in units of the electron rest mass energy
  //G4cout << "maxKE: " << maxKE << G4endl;
  return maxKE;
}

G4double SuNBetaDecay::GetRandomElectronKineticEnergy(const G4double maxKE){
  G4double REKE = maxKE * G4UniformRand(); // Calculates a random number between 0 and the maximum kinetic energy available to the electron in units of the electron rest mass energy
  return REKE;
}

G4double SuNBetaDecay::GetPhaseSpaceElectronKineticEnergyDistribution(const G4double KE, const G4double maxKE){
  G4double PSEKED = sqrt(KE * KE + 2.0 * KE * 1.0) * pow(maxKE - KE, 2.0) * (KE + 1.0); // Equation for the electron kinetic energy from Fermi's theory of beta decay in units of the electron rest mass energy. For example, see http://hyperphysics.phy-astr.gsu.edu/hbase/nuclear/beta2.html
  return PSEKED;
}

G4double SuNBetaDecay::GetPhaseSpaceElectronKineticEnergy(const G4double maxKE){

  G4double precision = 1000.0; // Determines how precisely maxValue is calculated below.

  std::vector<double> xValues(precision);
  std::vector<double> yValues(precision);
  
  for (int i=0; i<xValues.size(); i++){
    xValues[i] = (i/precision) * maxKE;
    yValues[i] = GetPhaseSpaceElectronKineticEnergyDistribution(xValues[i], maxKE); 
  }

  G4double maxValue = 0.0;
  for (int i=0; i<yValues.size(); i++){
    if (yValues[i] > maxValue){
      maxValue = yValues[i]; // Get the maximum value of the electron kinetic energy distribution
    }
  }
  
  G4double randomNumA, randomNumB, randomNumC;
  while (true){
    randomNumA = maxValue * G4UniformRand(); // Generate a random number between 0 and the maximum value of the electron kinetic energy distribution
    randomNumB = maxKE * G4UniformRand(); // Generate a random number between 0 and the maximum kinetic energy available to the electron in units of the electron rest mass energy
    randomNumC = GetPhaseSpaceElectronKineticEnergyDistribution(randomNumB, maxKE); // Run randomNumB through the electron kinetic energy distribution 
    if (randomNumA < randomNumC){
      G4double electronKE_mc2 = randomNumB; // The electron kinetic energy in units of the electron rest mass energy
      return electronKE_mc2; 
    }
  }
}

G4double SuNBetaDecay::GetApproximateFermiFunctionA(const G4double KE){
  G4double a, b, m, K, a0, b0, A, B, AFFA;
  if (zDN <= 56){
    a = 5.5465 * pow(10.0,-3.0);                           // ASR: pg 363; RDG: pg 10, eqn 3
    b = 76.929 * pow(10.0,-3.0);                           // ASR: pg 363; RDG: pg 10, eqn 3
  }
  if (zDN > 56){
    a = 1.2277 * pow(10.0,-3.0);                           // ASR: pg 363; RDG: pg 10, eqn 3
    b = 101.22 * pow(10.0,-3.0);                           // ASR: pg 363; RDG: pg 10, eqn 3
  }
  if (zDN < 16){
    m = 7.30 * pow(10.0,-2.0);
    K = 9.40 * pow(10.0,-1.0);
    A = m * zDN + K;                                       // ASR: pg 363; RDG: pg 10, eqn 2b
  }
  if (zDN >= 16){
    a0 = 404.56 * pow(10.0,-3.0);                          // ASR: pg 363;         RDG: pg 10, eqn 3
    b0 = 73.184 * pow(10.0,-3.0);                          // ASR: pg 363;         RDG: pg 10, eqn 3
    A = 1.0 + a0 * exp(b0*zDN);                            // ASR: pg 361, eqn 9a; RDG: pg 10, eqn 2a
  }
  B = a * zDN * exp(b*zDN);                                // ASR: pg 361, eqn 9b; RDG: pg 10, eqn 3
  AFFA = sqrt(A + (B / KE));                               // ASR: pg 361, eqn 8;  RDG: pg 10, eqn 1   // In the journal articles referenced above, the Fermi function is a function of the total energy. Here, it is a function of the kinetic energy: KE = (total energy) - (rest mass energy) = (total energy) - 1.0
  return AFFA;
}

G4double SuNBetaDecay::GetApproximateFermiFunctionB(const G4double KE){ // This code is adapted from the various websites listed above
  G4double e, p, u, s, y, a1, a2, AFFB;
  e = KE + 1.0;                                            // E = KE + mc^2
  p = sqrt(e * e - 1.0);                                   // E^2 = (pc)^2 + (mc^2)^2
  u = zDN / 137.0;                                         // The atomic number of the daughter nucleus multiplied by the fine structure constant
  s = sqrt(1.0 - u * u) - 1.0;                             // This is the S that is used in ASR: pg 360, eqn 7 and defined directly below that equation
  y = 2.0 * TMath::Pi() * u * e / p;                       // This is 2*pi*eta used in ASR: pg 360, eqn 3 and 7. Eta is defined on page 359
  a1 = u * u * e * e + p * p / 4.0;
  a2 = y / (1.0 - exp(-y));                                // This is F as used in ASR: pg 360, eqn 3.                                                                                              ==!!== There is an absolute value around this entire quantity in the websites I found, but I think it is is unneccessary ==!!==
  AFFB = pow(a1,s) * a2;                                   // This is F as used in ASR: pg 360, eqn 7 up to a constant (which is unimportant because it does not affect the shape of the function)
  return AFFB;
}

G4double SuNBetaDecay::GetApproximateFermiFunctionC(const G4double KE){ // Phys Rev C 81, 064314 (2010) Equation 12
  double e, p, u, s, y, R, AFFC;
  e = KE + 1.0;
  p = sqrt(e * e - 1.0);
  u = zDN / 137.0;
  s = sqrt(1.0 - u * u); // No -1.0 like there is in ApproximateFermiFunctionB
  y = u * e / p;
  R = (1.5*pow(10.0,-15.0))*pow(aDN,1.0/3.0) / (3.86159268*pow(10.0,-13.0)); // Need R to be dimensionless so divide by (hbar/(mc)) NOTE: R and (hbar/(mc)) must have same units
  AFFC = 2.0*(1.0+s) * pow(2.0*p*R,-2.0*(1.0-s)) * exp(TMath::Pi()*y) * std::norm(ComplexGammaFunction(s,y)) / std::norm(ComplexGammaFunction(2.0*s+1.0,0.0));
  return AFFC;
}

G4double SuNBetaDecay::GetCorrectedElectronKineticEnergy(const char flag, const G4double maxKE){

  char A = 'A'; // The possible flags for using which approximate Fermi function (these flag are used in PrimaryGeneratorAction.cc)
  char B = 'B';
  char C = 'C';

  G4double precision = 1000.0; // Determines how precisely maxValue is calculated below.

  std::vector<double> xValues(precision);
  std::vector<double> yValues(precision);

  for (int i=0; i<xValues.size(); i++){
    xValues[i] = (i/precision) * maxKE;
    if (flag == A){
      yValues[i] = GetApproximateFermiFunctionA(xValues[i]) * GetPhaseSpaceElectronKineticEnergyDistribution(xValues[i], maxKE); //  ==!!== Make sure that the Fermi function used here is the same as the one used directly above in the if statement and below in randomNumC ==!!==
      //G4cout << "using flag A the first time" << G4endl;
    }
    if (flag == B){
      yValues[i] = GetApproximateFermiFunctionB(xValues[i]) * GetPhaseSpaceElectronKineticEnergyDistribution(xValues[i], maxKE); //  ==!!== Make sure that the Fermi function used here is the same as the one used directly above in the if statement and below in randomNumC ==!!==
      //G4cout << "using flag B the first time" << G4endl;
    }
    if (flag == C){
      yValues[i] = GetApproximateFermiFunctionC(xValues[i]) * GetPhaseSpaceElectronKineticEnergyDistribution(xValues[i], maxKE); //  ==!!== Make sure that the Fermi function used here is the same as the one used directly above in the if statement and below in randomNumC ==!!==
      //G4cout << "using flag C the first time" << G4endl;
    }
  }

  G4double maxValue = 0.0;
  for (int i=0; i<yValues.size(); i++){
    if (yValues[i] > maxValue){
      maxValue = yValues[i]; // Get the maximum value of the electron kinetic energy distribution corrected by the Fermi function
    }
  }

  // Use the statistical "Acceptance-Rejection Method" to determine the electron kinetic energy
  G4double randomNumA, randomNumB, randomNumC;
  while (true){
    randomNumA = maxValue * G4UniformRand(); // Generate a random number between 0 and the maximum value of the electron energy distribution corrected by the Fermi function
    randomNumB = maxKE * G4UniformRand(); // Generate a random number between 0 and the maximum kinetic energy available to the electron in units of the electron rest mass energy
    if (flag == A){
      randomNumC = GetApproximateFermiFunctionA(randomNumB) * GetPhaseSpaceElectronKineticEnergyDistribution(randomNumB, maxKE); // Run randomNumB through the corrected electron kinetic energy distribution corrected by the Fermi function. ==!!== Make sure the Fermi function used here is the same as the one used directly above in the if statement and above in yValues[i] ==!!==
      //G4cout << "using flag A the second time" << G4endl; 
    }
    if (flag == B){
      randomNumC = GetApproximateFermiFunctionB(randomNumB) * GetPhaseSpaceElectronKineticEnergyDistribution(randomNumB, maxKE); // Run randomNumB through the corrected electron kinetic energy distribution corrected by the Fermi function. ==!!== Make sure the Fermi function used here is the same as the one used directly above in the if statement and above in yValues[i] ==!!==
      //G4cout << "using flag B the second time" << G4endl;
    }
    if (flag == C){
      randomNumC = GetApproximateFermiFunctionC(randomNumB) * GetPhaseSpaceElectronKineticEnergyDistribution(randomNumB, maxKE); // Run randomNumB through the corrected electron kinetic energy distribution corrected by the Fermi function. ==!!== Make sure the Fermi function used here is the same as the one used directly above in the if statement and above in yValues[i] ==!!==
      //G4cout << "using flag C the second time" << G4endl;
    }
    if (randomNumA < randomNumC){
      G4double electronKE_mc2 = randomNumB; // The electron kinetic energy in units of the electron rest mass energy
      return electronKE_mc2; 
    }
  }
}


 
