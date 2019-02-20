//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.							 *
//********************************************************************

//  Created: William Tireman - February 2019

#ifndef Hcal_h
#define Hcal_h

#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4ClassicalRK4.hh"

#include "NpolDetectorFactory.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhysicalVolume;

class Hcal : public NpolDetectorFactory {
  
public:
  Hcal();
  ~Hcal();
  
  void ConstructFakeHcal();
  void ConstructHCALarray();
  void ImprintPlate(G4AssemblyVolume *plate, G4LogicalVolume *motherLV,G4ThreeVector Tm,G4RotationMatrix Rm); 
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *motherLV);

  static G4double NpolAng, hcalZpos;

  static const G4Colour ScintColor,SteelColor,CopperColor;
  
private: 
  G4LogicalVolume *hcalLV;
  G4AssemblyVolume *hcalArray; 
};

#endif

