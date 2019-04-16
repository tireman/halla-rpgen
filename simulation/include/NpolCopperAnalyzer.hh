//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% NpolCopperAnalyzer.hh file %% //

// Created: William Tireman, November 2018

#ifndef Npol_CopperAnalyzer_h
#define Npol_CopperAnalyzer_h

#include "NpolDetectorFactory.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhysicalVolume;

class NpolCopperAnalyzer : public NpolDetectorFactory {
  
public:
  NpolCopperAnalyzer();
  ~NpolCopperAnalyzer();

  static G4double PosCopper, CopperThickness, CopperHeight, CopperWidth;

  void ConstructCopperAnalyzer();
  
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *motherLV);
  
private: 
  G4LogicalVolume *CopperAnalyzerLV;

};

#endif
