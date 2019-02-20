//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// && NpolScatteringChamber Header file %% //

// Created: Daniel Wilbern - March 2015

#ifndef Npol_ScatteringChamber_h
#define Npol_ScatteringChamber_h

#include "NpolDetectorFactory.hh"
#include "NpolTarget.hh"

class G4LogicalVolume;
class G4AssemblyVolume;
class G4VPhysicalVolume;

class NpolScatteringChamber : public NpolDetectorFactory {

public:
  NpolScatteringChamber();
  ~NpolScatteringChamber();

  void ConstructChamber();
  
  virtual G4String GetName();
  virtual void Place(G4LogicalVolume *motherLV);

  static G4double insideRadius;
  static G4double insideHeight;
  static G4double wallWidth;
  static G4double holeRadius;
  static G4double WindowThickness;
  static G4double shmsWindowAngle;
  static G4double shmsWindowDeltaAngle;
  static G4double npolWindowAngle;
  static G4double npolWindowDeltaAngle;
  static G4double npolWindowHeight;
  static G4double shmsWindowHeight;

private:
  
  G4LogicalVolume *innerChamberLV;
  G4LogicalVolume *chamberWallLV;
  NpolTarget *target;
};

#endif

