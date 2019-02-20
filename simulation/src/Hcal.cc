//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% Hcal.cc %%


#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "Hcal.hh"
#include "NpolPolarimeter.hh"

G4double Hcal::NpolAng = NpolPolarimeter::NpolAng;
G4double Hcal::hcalZpos = 830.0*cm; // 780*cm to front of HCAL; 50*cm to center

const G4Colour Hcal::ScintColor  = G4Colour(0.0,0.5,0.25);
const G4Colour Hcal::CopperColor  = G4Colour(0.0,1.0,0.0);
const G4Colour Hcal::SteelColor  = G4Colour(0.5,1.0,0.5);

Hcal::Hcal() {
  
  ConstructFakeHcal();
}

Hcal::~Hcal() {}

G4String Hcal::GetName() {
  return G4String("HCAL");
}

void Hcal::ConstructFakeHcal(){

  double hcalDetX = 15.*cm, hcalDetY = 15.*cm, hcalDetZ = 100.*cm;

  G4Box *hcalDet = new G4Box("hcalDet",hcalDetX/2,hcalDetY/2,hcalDetZ/2);

  G4LogicalVolume *hcalDetLV = new G4LogicalVolume(hcalDet,NpolMaterials::GetInstance()->GetMaterial("Scint"),"hcalDet",0,0,0);
  
  hcalDetLV->SetVisAttributes(ScintColor);
  hcalArray = new G4AssemblyVolume();
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  
  Ra.rotateX(0.0*deg); Ra.rotateY(0.0*deg); Ra.rotateZ(0.0*deg);
  
  // HCAL consists of 288 of these modules in a 12 by 24 array (wide x high)

  for(int i = 0; i < 24; i++){
	for(int j = 0; j < 12; j++){
	  double hcalXpos = 82.5*cm - j*15*cm;
	  double hcalYpos = 172.5*cm - i*15*cm;
	 
	  Ta.setX(hcalXpos); Ta.setY(hcalYpos); Ta.setZ(0.0*cm);
	  Tr = G4Transform3D(Ra,Ta);
	  hcalArray->AddPlacedVolume(hcalDetLV,Tr);

	}
  }
}

void Hcal::ConstructHCALarray(){
  //double hcalModuleX = 15.*cm, hcalModuleY = 15.*cm, hcalModuleZ = 100.*cm;
  /*double steelX = 7.25*cm, steelY = 15.*cm, steelZ = 1.25*cm;
  double scintX = 7.25*cm, scintY = 15.*cm, scintZ = 1.25*cm;
  double wlsX = 0.5*cm, wlsY = 15.*cm, wlsZ = 100.*cm;

  G4Box *steelPlate = new G4Box("steelPlate",steelX/2,steelY/2,steelZ/2);
  G4Box *scintPlate = new G4Box("scintPlate",scintX/2,scintY/2,scintZ/2);
  G4Box *wlsPlate = new G4Box("wlsPlate",wlsX/2,wlsY/2,wlsZ/2);

  G4LogicalVolume *steelPlateLV = new G4LogicalVolume(steelPlate,NpolMaterials::GetInstance()->GetMaterial("SSteel"),"steelPlateLV",0,0,0);
  G4LogicalVolume *scintPlateLV = new G4LogicalVolume(scintPlate,NpolMaterials::GetInstance()->GetMaterial("Scint"),"scintPlateLV",0,0,0);
  G4LogicalVolume *wlsPlateLV = new G4LogicalVolume(wlsPlate,NpolMaterials::GetInstance()->GetMaterial("Scint"),"wlsPlateLV",0,0,0);

  steelPlateLV->SetVisAttributes(SteelColor);
  scintPlateLV->SetVisAttributes(ScintColor);
  wlsPlateLV->SetVisAttributes(ScintColor);
  
  G4RotationMatrix ra, rm;
  G4ThreeVector ta, tm;
  
  G4AssemblyVolume *hcalModule = new G4AssemblyVolume();
  ta.setX(0.); ta.setY(0.); ta.setZ(0.);
  hcalModule->AddPlacedVolume(wlsPlateLV,ta, &ra);

  for(int i = 0; i < 40; i++){
	ta.setX(+(scintX + wlsX)/2); ta.setY(0.); ta.setZ(-(-1.25*cm + wlsZ/2 - scintZ*(2*i+1)/2 - steelZ*i));
	hcalModule->AddPlacedVolume(scintPlateLV,ta, &ra);
	ta.setX(+(steelX + wlsX)/2); ta.setY(0.); ta.setZ(-(-1.25*cm + wlsZ/2 - scintZ*(i+1) - steelZ*(2*i+1)/2));
	hcalModule->AddPlacedVolume(steelPlateLV,ta, &ra);
	ta.setX(-(scintX + wlsX)/2); ta.setY(0.); ta.setZ(-(wlsZ/2 - scintZ*(2*i+1)/2 - steelZ*i));
	hcalModule->AddPlacedVolume(scintPlateLV,ta, &ra);
	ta.setX(-(steelX + wlsX)/2); ta.setY(0.); ta.setZ(-(wlsZ/2 - scintZ*(i+1) - steelZ*(2*i+1)/2));
	hcalModule->AddPlacedVolume(steelPlateLV,ta, &ra);
  }
  // HCAL consists of 288 of these modules in a 12 by 24 array (wide x high)

  for(int i = 0; i < 24; i++){
	for(int j = 0; j < 12; j++){
	  double hcalXpos = 82.5*cm - j*15*cm;
	  double hcalYpos = 172.5*cm - i*15*cm;
	}
	}*/
  
}

void Hcal::Place(G4LogicalVolume *motherLV) {
  G4ThreeVector Ta;G4RotationMatrix Ra;G4Transform3D Tr;
  Ta.setX(-hcalZpos*sin(NpolAng*rad)); Ta.setY(0.0*cm); Ta.setZ(hcalZpos*cos(NpolAng*rad));
  Ra.rotateX(0.0*deg); Ra.rotateY(-NpolAng); Ra.rotateZ(0.0*deg);
  Tr = G4Transform3D(Ra,Ta);
  ImprintPlate(hcalArray,motherLV,Ta,Ra);
}


void Hcal::ImprintPlate(G4AssemblyVolume *plate, G4LogicalVolume *motherLV, G4ThreeVector Tm, G4RotationMatrix Rm) {
  
  G4Transform3D Tr;
  Tr = G4Transform3D(Rm,Tm);
  
  plate->MakeImprint(motherLV,Tr);
}
