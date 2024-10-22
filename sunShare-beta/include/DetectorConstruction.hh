
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"


extern  G4String detectorName[9];
extern  G4String detname[200];
extern int detno;
extern G4double X[10], Y[10], Z[10];

class G4VPhysicalVolume;
class G4Material;

class DetectorConstruction : public G4VUserDetectorConstruction
{
   public:
      DetectorConstruction();
      ~DetectorConstruction();

   public:
    G4VPhysicalVolume* Construct();

   private:
      void DefineMaterials();
      G4Material *NaI, *Al, *N78O21Ar1, *Cr20Ni8Fe76, *C2F4, *Si, *Cu3Zn2, *elCu, *C2H3Cl, *SiO2, *vacuum, *cardboard, *plastic;

      G4VPhysicalVolume* ConstructDetector();
};
 
#endif

