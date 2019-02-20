//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.							 *
//********************************************************************

// && NpolGEM Header file GEM detectors in NPOL

// Created: William Tireman - November 2018

#ifndef Npol_GEM_h
#define Npol_GEM_h

#include "globals.hh"
#include "G4ios.hh"
#include "G4Colour.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"

#include "NpolDetectorFactory.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhysicalVolume;

class NpolGEM : public NpolDetectorFactory {
  
public:
  NpolGEM();
  ~NpolGEM();

  
  //Logical volumes for Gem components
  G4LogicalVolume* GEMgapLV(){return GEMgap_log;};
  G4LogicalVolume* GEMKaptonLV(){return GEMKapton_log;};
  G4LogicalVolume* GEMFoilLV(){return GEMFoil_log;};
  G4LogicalVolume* GEMMylarLV(){return GEMMylar_log;};
  G4LogicalVolume* GEMdriftCathLV(){return GEMdriftCathode_log;};
  G4LogicalVolume* GEMReadoutLV(){return GEMReadout_log;};
  G4LogicalVolume* GEMdriftReadLV(){return GEMdriftReadout_log;}
  void BuildModule(G4double height, G4double width);
  void GEMLV();
  void ConstructDetector(G4AssemblyVolume* assemblyGEM,G4double numModules, G4double vertPosStart, G4double vertShift); 
  void ConstructINFNGEM();
  void ConstructUVaGEM();
  
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *motherLV);

private:
  G4double GEMFoil_x;
  G4double GEMFoil_y;
  G4double GEMFoil_z;
  
  G4double GEMKapton_x;
  G4double GEMKapton_y;
  G4double GEMKapton_z;
  
  G4double GEMgap_x;
  G4double GEMgap_y;
  G4double GEMgap_z;

  G4double GEMMylar_x;
  G4double GEMgapReadout_x;
  G4double z;

  G4double phi;
private:
  // Assembly Volume
  G4AssemblyVolume* assemblyINFNGEM;
  G4AssemblyVolume* assemblyUVaGEM;
  
  // Logical volumes
  G4LogicalVolume* GEMgap_log;
  G4LogicalVolume* GEMKapton_log; 
  G4LogicalVolume* GEMFoil_log;
  G4LogicalVolume* GEMMylar_log;
  G4LogicalVolume* GEMdriftCathode_log;
  G4LogicalVolume* GEMReadout_log;
  G4LogicalVolume* GEMdriftReadout_log;

  G4String GEMFoil_name = "GEMFoil";
  G4String GEMKapton_name = "GEMFoilKapton";
  G4String GEMgap_name = "GEMgap";
  G4String GEMReadout_name = "GEMreadout";
private:
  //VisAttributes of Gem components
  G4VisAttributes* GEMReadoutVisAtt; 
  G4VisAttributes* GEMdriftCathodeVisAtt;
  G4VisAttributes* GEMMylarVisAtt;
  G4VisAttributes* GEMFoilVisAtt;
  G4VisAttributes* GEMKaptonVisAtt;
  G4VisAttributes* GEMgapVisAtt; 
  G4VisAttributes* GEMdriftReadoutVisAtt;

  G4Box* GEMFoil;
  G4Box* GEMKapton;
  G4Box* GEMgap;
  G4Box* GEMMylar;
  G4Box* GEMdriftCathode;
  G4Box* GEMReadout;
  G4Box* GEMdriftR;

protected:
  static const G4Colour CopperColor; 
  static const G4Colour MylarColor;
  static const G4Colour NpolGasColor;
  static const G4Colour KaptonColor;

};

#endif
