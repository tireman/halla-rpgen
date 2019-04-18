//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	         *
//******************************************************************

// %% NpolCopperAnalyzer.cc %%

// Copper Analyzer construction file
// Created: William Tireman - November 2018

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolCopperAnalyzer.hh"
#include "NpolDipole2.hh"
#include "NpolPolarimeter.hh"

G4double NpolCopperAnalyzer::CopperThickness = 4.0*cm;  // thickness of the Cu analyzer
G4double NpolCopperAnalyzer::CopperHeight = 202*cm;
G4double NpolCopperAnalyzer::CopperWidth = 62.0*cm;
G4double NpolCopperAnalyzer::PosCopper = NpolPolarimeter::CuAnalyzerPos;

NpolCopperAnalyzer::NpolCopperAnalyzer() {
   ConstructCopperAnalyzer();
}

NpolCopperAnalyzer::~NpolCopperAnalyzer() {}

G4String NpolCopperAnalyzer::GetName() {
  return G4String("Copper Curtain");
}

// Construct a copper shield for in front of the polarimeter
void NpolCopperAnalyzer::ConstructCopperAnalyzer(){

  G4Box *CopperAnalyzer = new G4Box("CopperAnalyzer",CopperWidth/2,CopperHeight/2,CopperThickness/2);
  CopperAnalyzerLV = new G4LogicalVolume(CopperAnalyzer,NpolMaterials::GetInstance()->GetMaterial("Cu"),"CopperAnalyzerLV",0,0,0);
  G4VisAttributes *CopperAnalyzerVisAtt = new G4VisAttributes(G4Colour(1.0, 0.7, 0.2));
  CopperAnalyzerLV->SetVisAttributes(CopperAnalyzerVisAtt);
}


void NpolCopperAnalyzer::Place(G4LogicalVolume *motherLV) {
  G4double Angle = NpolPolarimeter::NpolAng;
  PlaceRectangular(CopperAnalyzerLV, motherLV, "CopperAnalyzer", -PosCopper*sin(Angle), 0, PosCopper*cos(Angle), 0.0*deg, -Angle, 0.0*deg);
 
}
