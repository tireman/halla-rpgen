//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.							 *
//********************************************************************

// Header file for the Polarimeter constructor

// Daniel Wilbern December/January 2014/15

#ifndef Npol_Polarimeter_h
#define Npol_Polarimeter_h

#include "G4SystemOfUnits.hh"

#include "NpolDetectorFactory.hh"
#include "NpolCopperAnalyzer.hh"
#include "NpolGEM.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhsicalVolume;

class NpolPolarimeter : public NpolDetectorFactory {
  
public:
  NpolPolarimeter();
  ~NpolPolarimeter();

  static G4double NpolAng;
  static G4double AnalyzerX, AnalyzerY, AnalyzerZ;
  static G4double LeftHodoscopeX,RightHodoscopeX, HodoscopeY, HodoscopeZ;
  static G4double HodoYPos, LeftHodoXPos, RightHodoXPos, HodoZPos;
  static G4double CuAnalyzerPos, CHAnalyzerPos;
  static const G4Colour CopperColor, SteelColor, ScintColor;
  
  void ConstructAnalyzerArray(G4LogicalVolume *motherLV);
  void ConstructHodoscopeArray(G4LogicalVolume *motherLV);
  void ConstructFakeGEM(G4LogicalVolume *motherLV);
  void ConstructPolarimeterFluxTagger(G4LogicalVolume *motherLV);
  
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *motherLV);
  void TranslateRotateAndPlace(G4LogicalVolume *polarimeterLV,G4LogicalVolume *motherLV,G4double rho,G4double phi,G4double z);
  
  G4AssemblyVolume *MakePlate(G4LogicalVolume *detLV,G4int numDets,G4double TmX,G4double TmY,G4double TmZ,G4double TmdX, G4double TmdY, G4double TmdZ);
  
  void ImprintPlate(G4AssemblyVolume *plate,G4LogicalVolume *motherLV,G4double TmX,G4double TmY,G4double TmZ,G4double RmZ);
  void ImprintPlate(G4AssemblyVolume *plate, G4LogicalVolume *motherLV,G4ThreeVector Tm,G4RotationMatrix Rm);
  
private:
  NpolCopperAnalyzer *copperAnalyzer;
  NpolGEM *gemDetectors;
  
};

#endif

