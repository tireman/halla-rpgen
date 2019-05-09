#ifndef __NpolSBSTargetBuilder_hh
#define __NpolSBSTargetBuilder_hh

#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

#include "NpolDetectorFactory.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhsicalVolume;
class G4DetectorConstruction;

class NpolSBSTargetBuilder: public NpolDetectorFactory { 
public:

  NpolSBSTargetBuilder();
  ~NpolSBSTargetBuilder();

  void BuildComponent(G4LogicalVolume *);
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *);

  void BuildCryoTarget(G4LogicalVolume *);
  void BuildGasTarget(G4LogicalVolume *);
   
  // EFuchey: 2017/02/10:  This function is now meant to build the cryotarget and target cell only.
  // This function takes as input the mother logical volume, a rotation matrix, and a 3-vector offset.
  void BuildStandardCryoTarget(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);

  // EFuchey: 2017/02/10: Added those functions to build scattering chamber separately from target,
  // and avoid, if possible, duplicates of the code actually building the target.
  void BuildStandardScatCham(G4LogicalVolume *);
    
  void SetTarget(G4String t){fTargType = t;}
  void SetTargLen(G4double len){ fTargLen = len;}
  void SetTargDen(G4double den){ fTargDen = den;} //Currently, fTargDen has NO effect!
  void SetTargDiameter(G4double D){ fTargDiameter = D; }
  void SetSchamFlag(int flag){ fSchamFlag = flag; }

  int GetSchamFlag() const { return fSchamFlag; }
  G4double GetTargLen() const { return fTargLen; }
  G4double GetTargDiameter() const { return fTargDiameter; }
 
  G4bool GetFlux() const { return fFlux; }
  void SetFlux(G4bool b){fFlux = b;}
  
private:
  G4double fTargLen;
  G4double fTargDen;
  G4double fTargDiameter; //diameter of cryotarget
  G4ThreeVector fTargPos; //Note: fTargPos currently has no effect!
  G4ThreeVector fTargDir; //Note: fTargDir currently has no effect!
  int fSchamFlag;

  G4bool fFlux;
  
  G4String fTargType,fExpType;

};

#endif//__NpolSBSTargetBuilder_hh
