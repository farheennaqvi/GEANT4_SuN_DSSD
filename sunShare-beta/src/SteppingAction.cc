////////////////////////////////////////////////////////////////////////////////////////////
// Description: This file tells the program what to do during each discrete step in       //
//   the simulation. For now the program simply reads in the energy deposited into each   //
//   detector and adds it up.                                                             //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////
 
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include <stdio.h>
#include <cstring>

SteppingAction::SteppingAction(DetectorConstruction* myDC):myDetector(myDC)
{ }


void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4VPhysicalVolume* currentVolume1 = aStep->GetPreStepPoint()-> GetPhysicalVolume();

  if (currentVolume1 != NULL)
  {  
    for(int i=0; i<nDetectors; i++)
    {
      if (currentVolume1->GetName() == detectorName[i])   
        energy_MeV[i] += aStep->GetTotalEnergyDeposit();
    }
  }
}
