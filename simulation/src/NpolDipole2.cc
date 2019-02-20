//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% NpolDipole2.cc %%

// Dipole 2 'BNL 48D48' Constructor file.  
// Created: William Tireman - November 2018
// Heavily modified from CGEN code

#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolDipole2.hh"
#include "NpolPolarimeter.hh"

// Field multiplier; if 2 mags, use value corresponding to desired integrated
// field strength; example: for 1 Tm use FM = 1.0 and for 4.3 Tm use FM = 4.3
// If 1 mag, use 2.0 for integrated 1.0 Tm or 4.0 for 2.0 Tm (double FM for 
// same effect) Note: Generally we won't use Dipole 2 by itself
G4double NpolDipole2::inch = 2.54*cm;
G4double NpolDipole2::FM = 2.0;    //Integrate field length  (Tesla * meters) 
G4double NpolDipole2::dipole2FieldY = FM*0.819672*tesla;   // Inverted 1.22 meters multiplier

G4double NpolDipole2::NpolAng = NpolPolarimeter::NpolAng;
G4double NpolDipole2::yokeLength = 48.0*inch; // 1.22*m;
G4double NpolDipole2::gapWidth = 48.0*inch; //1.22*m;
G4double NpolDipole2::gapLength = 48.0*inch; //1.22*m;
G4double NpolDipole2::gapHeight = 18.5*inch; //0.4699*m;
G4double NpolDipole2::PosD2 = 292.0*cm; // RP-GEN version
G4double NpolDipole2::fieldClampHeight = 280.82*cm;
G4double NpolDipole2::fieldClampWidth = 160.02*cm;
G4double NpolDipole2::fieldClampThick = 10.16*cm;
G4double NpolDipole2::fieldClampInheight = 70.0*cm;
G4double NpolDipole2::fieldClampInwidth = 91.10*cm;

G4double NpolDipole2::ClampOffSetF = 0.5*gapLength+0.5*fieldClampThick+34.04*cm;
G4double NpolDipole2::ClampOffSetR = 99.70*cm;

NpolDipole2::NpolDipole2() {
 
}

NpolDipole2::~NpolDipole2() {}

G4String NpolDipole2::GetName() {
  return G4String("Dipole 2");
}

void NpolDipole2::ConstructDipole2Yoke(G4LogicalVolume *motherLV) {

  G4double HOIoutsideY = 371.43*cm; G4double HOIoutsideX = 231.29*cm; G4double HOIoutsideZ = 121.*cm;
  G4double HOIinsideY = 187.0*cm; G4double HOIinsideX = 47.5*cm; G4double HOIinsideZ = 121.*cm + 1*cm;
  G4Box *hunkOfIron = new G4Box("hunkOfIron",HOIoutsideX/2, HOIoutsideY/2, HOIoutsideZ/2);
  G4Box *insideCutout =
	new G4Box("insideCutout",HOIinsideX/2, HOIinsideY/2, HOIinsideZ/2);
  G4Box *beamlineCutout =
	new G4Box("beamlineCutout", (46.38/2)*cm, (33.2/2)*cm, (123/2)*cm);
  G4Trap *beamlineTrap =
	new G4Trap("beamlineTrap", 36.57*cm, 22.21*cm, 25.72*cm/2, 33.2*cm/2, 122*cm/2);
  
  G4SubtractionSolid *step1 =
	new G4SubtractionSolid("step1",hunkOfIron,insideCutout);
  G4SubtractionSolid *step2 =
	new G4SubtractionSolid("step2",step1,beamlineCutout,0,G4ThreeVector(92.465*cm,0, 0));
 
  G4SubtractionSolid *Dipole2Yoke =
	new G4SubtractionSolid("Dipole2Yoke",step2,beamlineTrap, 0, G4ThreeVector((69.285*cm),0,0));
  
  Dipole2YokeLV = 
	new G4LogicalVolume(Dipole2Yoke, NpolMaterials::GetInstance()->GetMaterial("Fe"),"Dipole2YokeLV",0,0,0);
  
  G4VisAttributes *Dipole2YokeVisAtt=  new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  Dipole2YokeLV->SetVisAttributes(Dipole2YokeVisAtt);
  
   PlaceRectangular(Dipole2YokeLV, motherLV, "Dipole2", -(PosD2*sin(NpolAng)), 0.0*cm, +(PosD2*cos(NpolAng)), 0*deg, -NpolAng, 0*deg);
}

void NpolDipole2::ConstructDipole2LeadInsert(G4LogicalVolume *motherLV){
  
  G4Box *leadBox =
	new G4Box("leadBox",  (68.0/2)*cm, (32.0/2)*cm, (60/2)*cm);
  leadInsertLV =
	new G4LogicalVolume(leadBox,NpolMaterials::GetInstance()->GetMaterial("Pb"),"leadInsertLV",0,0,0);
  
  G4VisAttributes *leadVisAtt=  new G4VisAttributes(G4Colour(1.0,0.5,1.0));
  leadInsertLV->SetVisAttributes(leadVisAtt);
  G4double xOffset = 81.0*cm; G4double zOffset = 30.0*cm;
  G4double xPos
	= xOffset*cos(-NpolAng) + zOffset*sin(-NpolAng) - (PosD2)*sin(NpolAng);
  G4double zPos
	= -xOffset*sin(-NpolAng) + zOffset*cos(-NpolAng)+ (PosD2)*cos(NpolAng);
 
  PlaceRectangular(leadInsertLV, motherLV, "LeadInsert", xPos, 0.0*cm, zPos,0*deg, -NpolAng, 0*deg);

}

void NpolDipole2::ConstructDipole2CoilPack(G4LogicalVolume *motherLV){

  G4double xPos, zPos;
  G4double xOffset = -(12.5*cm+(47.5*cm)/2); G4double zOffset = 0.0*cm;
  G4Box *hunkOfCopper = new G4Box("hunkOfCopper", 71.8*cm/2, 187*cm/2, 165*cm/2);
  G4Box *horCenterCut = new G4Box("horCenterCut", (41.0*cm/2), (127.5*cm/2), 166.1*cm/2);
  G4Box *transCenterCut = new G4Box("transCenterCut", 50.0*cm/2, 188.0*cm/2, 123.6*cm/2);

  G4SubtractionSolid *step1 =
	new G4SubtractionSolid("step1",hunkOfCopper,horCenterCut, 0, G4ThreeVector(15.5*cm, 0*cm, 0*cm));
  G4SubtractionSolid *Dipole2SaddleCoil =
	new G4SubtractionSolid("Dipole2SaddleCoil",step1,transCenterCut, 0, G4ThreeVector(-12.5*cm, 0*cm, 0*cm));
															   
  Dipole2SaddleCoilLV =
	new G4LogicalVolume(Dipole2SaddleCoil,NpolMaterials::GetInstance()->GetMaterial("Cu"),"Dipole2SaddleCoilLV", 0,0,0);
  G4VisAttributes *CuMat = new G4VisAttributes(G4Colour(0.50,0.50,0.05));
  Dipole2SaddleCoilLV->SetVisAttributes(CuMat);

  xPos = xOffset*cos(-NpolAng) + zOffset*sin(-NpolAng)- (PosD2)*sin(NpolAng);
  zPos= -xOffset*sin(-NpolAng) + zOffset*cos(-NpolAng)+ (PosD2)*cos(NpolAng);

  PlaceRectangular(Dipole2SaddleCoilLV,motherLV,"Dipole2SaddleCoil", xPos,0*cm,zPos,0*deg,-NpolAng,0*deg);

  G4Box *hunkOfCopper2 = new G4Box("hunkOfCopper2", 139.8*cm/2, 187*cm/2, 165*cm/2);
  G4Box *horCenterCut2 = new G4Box("horCenterCut2", (140.5*cm/2), (143.0*cm/2), 166.0*cm/2);
  G4Box *transCenterCut2 = new G4Box("transCenterCut2", 96.05*cm/2, 187.0*cm/2, 124.0*cm/2);
  
  G4SubtractionSolid *step2 =
	new G4SubtractionSolid("step2",hunkOfCopper2,horCenterCut2, 0, G4ThreeVector(0, 0, 0));
  G4SubtractionSolid *Dipole2RaceCoil =
	new G4SubtractionSolid("Dipole2RaceCoil",step2,transCenterCut2, 0, G4ThreeVector(0, 0, 0));
  
  Dipole2RaceCoilLV =
	new G4LogicalVolume(Dipole2RaceCoil,NpolMaterials::GetInstance()->GetMaterial("Cu"),"Dipole2RaceCoilLV", 0,0,0);
  Dipole2RaceCoilLV->SetVisAttributes(CuMat);
  
  xOffset = (231.29*cm+47.0*cm)/4;
  xPos = xOffset*cos(-NpolAng) + zOffset*sin(-NpolAng)- (PosD2)*sin(NpolAng);
  zPos= -xOffset*sin(-NpolAng) + zOffset*cos(-NpolAng)+ (PosD2)*cos(NpolAng);
  PlaceRectangular(Dipole2RaceCoilLV, motherLV, "Dipole2RaceCoil", xPos,0.0*cm, zPos, 0*deg, -NpolAng, 0*deg);
  
}

void NpolDipole2::ConstructDipole2FieldClamp(G4LogicalVolume *motherLV){
  G4double xPos, zPos;
  G4double xOffset = 0.0*cm; G4double zOffset = -ClampOffSetF;
  
  G4Box *Slab = new G4Box("Slab", fieldClampWidth/2, fieldClampHeight/2, fieldClampThick/2);
  G4Box *beamlineCutout = new G4Box("beamlineCutout", (fieldClampInwidth)/2, (fieldClampInheight)/2, (fieldClampThick+0.001*m)/2);

  double cutOffset = 0.5*(fieldClampWidth-fieldClampInwidth);
  G4SubtractionSolid *frontFieldClamp =
	new G4SubtractionSolid("frontFieldClamp",Slab,beamlineCutout,0,G4ThreeVector(cutOffset,0,0));
  
  frontFieldClampLV = 
	new G4LogicalVolume(frontFieldClamp, NpolMaterials::GetInstance()->GetMaterial("Fe"),"frontFieldClampLV", 0,0,0);
  G4VisAttributes *Clamp = new G4VisAttributes(G4Colour(0.1,0.5,0.1));
  frontFieldClampLV->SetVisAttributes(Clamp);
  
  xPos = xOffset*cos(-NpolAng) + zOffset*sin(-NpolAng)- (PosD2)*sin(NpolAng);
  zPos= -xOffset*sin(-NpolAng) + zOffset*cos(-NpolAng)+ (PosD2)*cos(NpolAng);
  PlaceRectangular(frontFieldClampLV,motherLV,"frontFieldClamp",xPos,0.0*cm,zPos,0.0*deg,-NpolAng,0.0*deg);


  /*G4Box *SlabR = new G4Box("SlabR", 292*cm/2, 267*cm/2, 20*cm/2);
  G4Box *centerHoleR = new G4Box("centerHoleR", 130*cm/2, 46*cm/2, 21*cm/2);
  G4Box *beamlineCutoutR = new G4Box("beamlineCutoutR", 40*cm/2, 50.0*cm/2, 21*cm/2);
  G4Box *chopSlab = new G4Box("chopSlab", 180*cm/2, 268*cm/2, 10.1*cm/2);
  
  G4SubtractionSolid *step1R = new G4SubtractionSolid("step1R", SlabR, centerHoleR, 0, G4ThreeVector(0,-17.5*cm,0));
  G4SubtractionSolid *step2R = new G4SubtractionSolid("step2R", step1R, chopSlab, 0 , G4ThreeVector(0, 0, +5.05*cm));
  G4SubtractionSolid *rearFieldClamp =
	new G4SubtractionSolid("rearFieldClamp",step2R,beamlineCutoutR,0,G4ThreeVector(0,82.5*cm,0));

  rearFieldClampLV = new G4LogicalVolume(rearFieldClamp, NpolMaterials::GetInstance()->GetMaterial("Fe"),"rearFieldClampLV", 0,0,0);
  rearFieldClampLV->SetVisAttributes(Clamp);
  
  PlaceRectangular(rearFieldClampLV,motherLV,"rearFieldClamp",+17.5*cm,0.0*cm,+ClampOffSetR,0.0*deg,0.0*deg,-90.0*deg);*/
  
}

void NpolDipole2::ConstructDipole2Field(G4LogicalVolume *motherLV){
  G4double xPos, zPos;
  G4double xOffset = 0.0*cm; G4double zOffset = 0.0*cm;
  
  // Generate the magnetic field volume
  G4Box *bigFieldBox = new G4Box("bigFieldBox",47.0*cm/2, 143*cm/2, 122*cm/2);
  G4Box *fieldSub =
	new G4Box("fieldSub",23.6*cm/2, 11.50*cm/2, 123*cm/2);  // tweaked to get to work.

  G4SubtractionSolid *fieldStep1 = new G4SubtractionSolid("fieldStep1",bigFieldBox,fieldSub,0, G4ThreeVector(-11.75*cm, -(143*cm/2 - 2.5*cm), 0));
  G4SubtractionSolid *Dipole2Field = new G4SubtractionSolid("Dipole2Field", fieldStep1,fieldSub,0, G4ThreeVector(-11.75*cm, (143*cm/2 - 2.5*cm), 0));
  
  Dipole2FieldLV = 
	new G4LogicalVolume(Dipole2Field, NpolMaterials::GetInstance()->GetMaterial("Vacuum"),"Dipole2FieldLV", 0,0,0);
  G4VisAttributes *Field = new G4VisAttributes(G4Colour(0.5,0.7,0.2));
  Dipole2FieldLV->SetVisAttributes(Field);
  Dipole2FieldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  xPos = xOffset*cos(-NpolAng) + zOffset*sin(-NpolAng)- (PosD2)*sin(NpolAng);
  zPos= -xOffset*sin(-NpolAng) + zOffset*cos(-NpolAng)+ (PosD2)*cos(NpolAng);
  PlaceRectangular(Dipole2FieldLV, motherLV, "Field", xPos, 0.0*cm, zPos, 0.0*deg,-NpolAng,0.0*deg);
  
}

void NpolDipole2::Place(G4LogicalVolume *motherLV) {
   ConstructDipole2Yoke(motherLV);
   ConstructDipole2LeadInsert(motherLV);
   ConstructDipole2CoilPack(motherLV);
   ConstructDipole2FieldClamp(motherLV);
   ConstructDipole2Field(motherLV);
  // End of Dipole 2 (BNL 48D48) construction.  May it rest in place.
}

void NpolDipole2::ConstructSDandField() {
  // Generate the Magnetic Field for BNL
  G4TransportationManager* tmanMagField = G4TransportationManager::GetTransportationManager();
  tmanMagField -> GetPropagatorInField() -> SetLargestAcceptableStep(1*mm);
  
  G4UniformMagField *magField = new G4UniformMagField(G4ThreeVector(dipole2FieldY, 0.0, 0.0));
  G4Mag_EqRhs *fEqMagField = new G4Mag_UsualEqRhs(magField);
  G4double minStepMagneticField = 0.0025*mm;
  G4FieldManager *fieldManMagField = new G4FieldManager(magField);
  
  G4MagIntegratorStepper *stepperMagField = new G4ClassicalRK4(fEqMagField);
  fieldManMagField -> SetDetectorField(magField);
  
  G4ChordFinder *fChordFinder = new G4ChordFinder(magField, minStepMagneticField, stepperMagField);
  fieldManMagField->SetChordFinder(fChordFinder);
  
  Dipole2FieldLV->SetFieldManager(fieldManMagField, true);
}

