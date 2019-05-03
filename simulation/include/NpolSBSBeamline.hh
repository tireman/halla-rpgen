#ifndef __NpolSBSBeamline_h
#define __NpolSBSBeamline_h

#include "G4SystemOfUnits.hh"
#include "NpolDetectorFactory.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhsicalVolume;

class NpolSBSBeamline : public NpolDetectorFactory {

public:
  NpolSBSBeamline();
  ~NpolSBSBeamline();
  
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *motherLV);

  static G4int ScattChamberFlag, LeadOption,fBeamlineConf;
  static G4String TargType, ExpType;
  static G4double beamheight, ent_len, swallrad_inner, swallrad;
  static G4double ent_rin, ent_rou,gapheight,shieldblock3_height;
  static G4double f48D48depth, f48D48width, f48D48height;
  static G4double f48D48ang, f48D48dist, f48D48_fieldclamp_config;
  static G4double leadstart,leadend,magleadlen;
private:

  void ConstructGEpLead();
  void ConstructGMnBeamline(G4LogicalVolume *);
  void ConstructGMnLead();
  void ConstructCommonExitBeamline(G4LogicalVolume *);

private:

  G4LogicalVolume *entLog, *entvacLog, *entLog_cut, *entvacLog_cut;
  G4LogicalVolume *extLog, *extvacLog, *floorLog;
  G4LogicalVolume *clamp_cone1_log,*clamp_cone2_log,*clamp_cone3_log;
  G4LogicalVolume *leadinmag_log,*windowshield_log,*shieldblock3_log;
  G4LogicalVolume *shield_cone1_log,*leadcone1_log,*leadshield2_log;
  G4LogicalVolume *lead_wall1_log,*lead_wall2_log;
  
  G4double fHCALdist  = 17.0*m;
  G4double fHCALvertical_offset = 0.0*cm;
  G4double fRICHdist  = 15.0*m;
  

};

#endif//__NpolSBSBeamline_h
