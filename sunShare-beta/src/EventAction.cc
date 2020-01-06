////////////////////////////////////////////////////////////////////////////////////////////
// Description: This file tells the program what to do before and after each              //
//   gamma-ray cascade. For now the program does the following:                           //
//                                                                                        //
//   BEFORE: - Calculate the gamma-ray cascade to run.                                    //
//           - Set energies and multiplicities back to zero.                              //
//                                                                                        //
//   AFTER:  - Calculate the energy deposited based on your detector resolution function. //
//           - Calculate the multiplicity of detectors hit.                               //
//           - Fill the ROOT Tree.                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <stdio.h>
#include "Randomize.hh"

EventAction::EventAction()
{}

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  //status bar that prints out run number to screen
  if ((evt->GetEventID()+1) % 10000 == 0)                 
   G4cout << ">>> Event " << evt->GetEventID()+1 << G4endl;

 // figures out what gamma-ray cascade to run
  int j = evt->GetEventID()+1; 
  for (int i=0; i<nCascades; i++) 
  {
     if (i==0 && j<inputCutoff[0])
       cascade=i;
     else if (i>0 && j<inputCutoff[i] && j>=inputCutoff[i-1])
       cascade=i;
  }

  //set everything back to zero
  for(int i=0; i<nDetectors; i++) 
   { 
    energy_MeV[i] = 0.0*eV;
    energy_keV[i] = 0.0*keV;
    sigma[i] = 0.0*keV;
    energy[i] = 0.0*keV;
   }
  energy_tot=0.0*keV;
  mult=0;
}


void EventAction::EndOfEventAction(const G4Event* evt)
{

  // calculate energy and multiplicity to save to ROOT file
  for(int i=0; i<8; i++)
   {
     double threshold = 40.0; //experimental threshold in keV

     if (energy_MeV[i]*1000.0 > threshold)
	   {	
	     energy_keV[i]=1000.0*energy_MeV[i];

       //detector resolution function
       sigma[i]=-5.59375e-15*pow(energy_keV[i],4.0)
                +1.85975e-10*pow(energy_keV[i],3.0)
                -2.47836e-6 *pow(energy_keV[i],2.0)
                +2.33408e-2 *energy_keV[i]
                +7.00328;

       energy[i] = G4RandGauss::shoot(energy_keV[i],sigma[i]);

       energy_tot += energy[i];
       ++mult;
     }
     else energy[i]=0.0; 
   }

  // calculate energy for plug to save in root
  double threshold_plug = 10.0; // change to exp threshold 
     if (energy_MeV[8]*1000.0 > threshold_plug)
	   {	
	     energy_keV[8]=1000.0*energy_MeV[8];

       //detector resolution function for NaI
       sigma[8]=0.009 *energy_keV[8];
               

       energy[8] = G4RandGauss::shoot(energy_keV[8],sigma[8]);

	   }
     else energy[8]=0.0; 




    Ts->Fill();
}
