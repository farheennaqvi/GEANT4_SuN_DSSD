//   
//

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"


extern int nEnergies, cascade;
extern float inputArray[5100][50];
extern double electronKE_mc2;
extern double electronKE_keV;

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event *anEvent);

  private:
    G4ParticleGun* particleGun;
    G4ParticleGun* particleGun2;

};

#endif

