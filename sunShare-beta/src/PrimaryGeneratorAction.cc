////////////////////////////////////////////////////////////////////////////////////////////
// Description: This file tells the program what particle to shoot, where to shoot it out //
//   from, what angle to shoot it at, and at what energy the particle has.                //
//                                                                                        //
// Most of the code added below was written to include the Doppler effect for particles   //
// emitted in-flight. The kinetic energy of the moving source is defined by the beamE     //
// variable. Most of the time beamE=0 and the particle is emitted at rest.                //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

#include "PrimaryGeneratorAction.hh"
#include "SuNBetaDecay.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh" 
#include "G4ParticleDefinition.hh"      
#include "G4UImanager.hh"

#include "globals.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <iomanip>

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;

  // Setting the particle to a gamma ray
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4String pName;
  G4ParticleDefinition *particle = pTable->FindParticle(pName="gamma");

  particleGun = new G4ParticleGun(n_particle);
  particleGun->SetParticleDefinition(particle);

  G4ParticleDefinition *particle2 = pTable->FindParticle(pName="e-");
  particleGun2 = new G4ParticleGun(n_particle);
  particleGun2->SetParticleDefinition(particle2);


}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete particleGun2;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)    
{
  //variables
  G4double x, y, z;                       // position of source
  G4double costheta, phi;                 // angle in CoM frame
  G4double cosLabTheta, sinLabTheta;      // angle in Lab frame
  G4double pxLab, pyLab, pzLab;           // momentum vector in Lab frame
  G4double dopplerEffect, gmma, velocity; // doppler variables
  G4double beamE, mamu;                   // energy variables

  mamu= 931.5 *MeV;
  beamE=0.0; //*MeV/u

  for(int i=1; i<=nEnergies+1; i++)
   {
   // Generate random angle to shoot particle at
    costheta=2.0*G4UniformRand()-1.0;         // generates a random value for costheta between -1 and 1
    phi= twopi*G4UniformRand();               // generates a random angle between 0 and 2pi

    // Now we start necessary calcs for Doppler effects (I believe this was initially developed by Shea Mosby)
     gmma = ((beamE*MeV)/mamu) + 1.;          // E_kinetic = gamma*mass*c*c - mass*c*c
     velocity = sqrt(1. - (1./(gmma*gmma)));  // gamma = 1/sqrt(1-(v/c)^2)

   // Move to the lab frame (we assume here that the beam axis is the z axis)
     cosLabTheta = (costheta + velocity)/(1. + velocity*costheta);   // relativistic effects
     sinLabTheta=std::sqrt(1.0-cosLabTheta*cosLabTheta);             // sin^2 + cos^2 = 1
     pxLab=sinLabTheta * std::cos(phi);                              // calculates x = sintheta * cosphi
     pyLab=sinLabTheta * std::sin(phi);                              // calculates y = sintheta * sinphi
     pzLab=cosLabTheta;                                              // calculates z = costheta
     dopplerEffect = 1./(gmma*(1. - velocity*cosLabTheta));          // f = 1/gamma * 1/(1-v/c*costheta) * f0

   // Sets the original postion of the particle
     x=0.0*cm;
     y=0.0*cm;
     z=0.0*cm;

   // Sets the particle's postition, momentum vector, and energy with the Doppler boost
     particleGun->SetParticlePosition(G4ThreeVector(x, y, z));

     particleGun->SetParticleMomentumDirection(G4ThreeVector(pxLab, pyLab, pzLab));

     if (i<=nEnergies){ // Set the energies for the gamma rays
       particleGun->SetParticleEnergy(dopplerEffect*inputArray[cascade][i]*keV);
       // To generate an event
       particleGun->GeneratePrimaryVertex(anEvent);
       /*
	 G4cout << "Cascade: " << cascade << "  " << "Particle Definition for particle " << i << " is: " << particleGun->GetParticleDefinition() << G4endl;
	 G4cout << "Cascade: " << cascade << "  " << "Particle Energy for particle " << i << " is: " << particleGun->GetParticleEnergy() << G4endl;
	 G4cout << "Cascade: " << cascade << "  " << "Particle Number for particle " << i << " is: " << particleGun->GetNumberOfParticles() << G4endl;
       */
     }
     
     if (i==nEnergies+1){  // Set the energy for the electron

       SuNBetaDecay Zr101BetaDecay;
       G4double maxKE = Zr101BetaDecay.GetMaximumElectronKineticEnergy();
       //electronKE_mc2 = 0.0;
       //electronKE_mc2 = Zr101BetaDecay.GetRandomElectronKineticEnergy(maxKE);
       //electronKE_mc2 = Zr101BetaDecay.GetPhaseSpaceElectronKineticEnergy(maxKE);
       electronKE_mc2 = Zr101BetaDecay.GetCorrectedElectronKineticEnergy('B',maxKE);
       electronKE_keV = electronKE_mc2 * (Zr101BetaDecay.electronRestMassEnergy);
       particleGun2->SetParticlePosition(G4ThreeVector(x, y, z));
       particleGun2->SetParticleMomentumDirection(G4ThreeVector(pxLab, pyLab, pzLab));
       particleGun2->SetParticleEnergy(electronKE_keV*keV);
       particleGun2->GeneratePrimaryVertex(anEvent);
       /*
	 G4cout << "electronKE_mc2: " << electronKE_mc2 << G4endl;
	 G4cout << "electronKE_keV: " << electronKE_keV << G4endl;
	 G4cout << "Cascade: " << cascade << "  " << "Particle Definition for particle " << i << " is: " << particleGun2->GetParticleDefinition() << G4endl;
	 G4cout << "Cascade: " << cascade << "  " << "Particle Energy for particle " << i << " is: " << particleGun2->GetParticleEnergy() << G4endl;
	 G4cout << "Cascade: " << cascade << "  " << "Particle Number for particle " << i << " is: " << particleGun2->GetNumberOfParticles() << G4endl;
	 G4cout << "************************************************************" << G4endl;
       */


   }

}
