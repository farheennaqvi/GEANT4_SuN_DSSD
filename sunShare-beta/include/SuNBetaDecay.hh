#include "globals.hh"

class SuNBetaDecay
{
public:
  G4int zDN; // The atomic number of the daughter nucleus in the beta decay (D = daughter, N = nucleus)

  G4int aDN; // The mass number of the daughter nucleus in the beta decay (D = daughter, N = nucleus)       

  G4double Qvalue; // The ground-state-to-ground-state Q value in keV for the beta decay. This value is from the NNDC.

  G4double electronRestMassEnergy; // The rest mass energy of the electron in keV.

  SuNBetaDecay()
  {
    zDN = 30;                           //Cu74
    aDN = 74;
    Qvalue = 9751.0;
    electronRestMassEnergy = 510.99891;
  }
  
  ~SuNBetaDecay()
  {}
  
  std::complex<G4double> ComplexGammaFunction(const G4double a, const G4double b);

  G4double GetMaximumElectronKineticEnergy();

  G4double GetRandomElectronKineticEnergy(const G4double KE);

  G4double GetPhaseSpaceElectronKineticEnergyDistribution(const G4double KE, const G4double maxKE);

  G4double GetPhaseSpaceElectronKineticEnergy(const G4double maxKE);
           
  G4double GetApproximateFermiFunctionA(const G4double KE);

  G4double GetApproximateFermiFunctionB(const G4double KE);

  G4double GetApproximateFermiFunctionC(const G4double KE);

  G4double GetCorrectedElectronKineticEnergy(const char flag, const G4double maxKE);

private:
};
