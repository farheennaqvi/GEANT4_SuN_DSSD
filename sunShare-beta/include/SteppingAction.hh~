// 
// 

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"


extern int nDetectors;
extern G4double energy_MeV[8]; 
extern G4String detectorName[8];


class DetectorConstruction;
class G4Track;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction* myDC);
    virtual ~SteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

  private:
    DetectorConstruction* myDetector;
};

#endif

