// 
// 

#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"


extern G4double energy_MeV[8], energy_keV[8], energy[8], sigma[8], energy_tot;
extern int mult, cascade, nCascades, nDetectors;
extern float inputCutoff[500];


class G4Event;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
};

#endif

