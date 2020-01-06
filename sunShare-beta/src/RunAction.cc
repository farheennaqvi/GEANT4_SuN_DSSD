///////////////////////////////////////////////////////////////////////////////////
// Description: This file tells the program what to do before and after running. //
//   For now the program does the following:                                     //
//                                                                               //
//   BEFORE: - Creates the ROOT tree and branches                                //
//                                                                               //
//   AFTER:  - Write the data to the ROOT tree                                   //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4String.hh"

#include <stdio.h>
#include <sstream>
#include <iostream>
#include "globals.hh"

RunAction::RunAction()
{
  runIDcounter = 0;
}


RunAction::~RunAction()
{}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  ((G4Run *)(aRun))->SetRunID(runIDcounter++);

  G4cout << "Run " << aRun->GetRunID() << " start." << G4endl;

  srand(time(0));

  //creates the file and tree for the .root output
  newfile = new TFile("rootfiles/out.root","recreate");
  Ts = new TTree("Ts","output from geant");

  G4String branchName;
  G4String branchNameD;

  for(int i=0; i<nDetectors; i++)
   {
     branchName="ene"+detectorName[i];
     branchNameD="ene"+detectorName[i]+"/D";
 
     ebranch = Ts->Branch(branchName,&energy[i],branchNameD);
   }

   ebranch = Ts->Branch("eneAll",&energy_tot,"eneAll/D");
   ebranch = Ts->Branch("multi",&mult,"multi/I");

  if(G4VVisManager::GetConcreteInstance())
   {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/clear/view");
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");
   }
}


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "Run " << aRun->GetRunID() << " ended." << G4endl;

  if(G4VVisManager::GetConcreteInstance())
   {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
   }

  Ts->Write(); //write the tree to file
}












