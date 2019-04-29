//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% NpolPolarimeter.cc %%

// Modified from the C-GEN code, this is for the Recoil Polarimetry Experiment
// in Hall A.

// Polarimeter construction file
// Created: William Tireman - November 2018
// This code currently (November 2018) generates NPOL volumes with
// names of the following scheme av_xxx_impr_yyy_VolName_pv_zzz
// xxx is the assembly volume number (runs from 1 to # of assemblies)
// yyy is the imprint number of the assembly for multiple placements
//     of the same assembly (runs from 1 to # of imprints)
// zzz is the physical volume number for detectors in an assembly
//     (runs from 0 to the # of detectors in the assembly)
// VolName is the physical volume name given when it is created

#include <string>
#include <cstring>
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolPolarimeter.hh"
#include "NpolGEM.hh"
#include "NpolAnalysisManager.hh"
#include "NpolDetectorFactory.hh"
#include "NpolDetectorConstruction.hh"
#include "NpolCopperAnalyzer.hh"

G4double NpolPolarimeter::NpolAng = 24.7*deg;

G4double NpolPolarimeter::AnalyzerX = 4.0*cm;
G4double NpolPolarimeter::AnalyzerY = 4.0*cm;
G4double NpolPolarimeter::AnalyzerZ = 25.*cm;
G4double NpolPolarimeter::vertAnalyzerX = 10.0*cm;
G4double NpolPolarimeter::vertAnalyzerY = 100.0*cm;
G4double NpolPolarimeter::vertAnalyzerZ = 10.0*cm;

G4double NpolPolarimeter::LeftHodoscopeX = 0.3*cm;
G4double NpolPolarimeter::RightHodoscopeX = 3.0*cm;
G4double NpolPolarimeter::HodoscopeY = 8.6*cm;
G4double NpolPolarimeter::HodoscopeZ = 50.0*cm;

// Cu Analyzer from target; also controls z-position of all polarimeter parts
G4double NpolPolarimeter::CuAnalyzerPos = 478.0*cm; // to the center of Cu Analyzer
G4double NpolPolarimeter::CHAnalyzerPos = CuAnalyzerPos + 63.82*cm; // offset from Cu
G4double NpolPolarimeter::LeftHodoXPos = 127.15*cm; // x-direction (left right of beam)
G4double NpolPolarimeter::RightHodoXPos = 128.35*cm; // x-direction (left right of beam)
G4double NpolPolarimeter::HodoYPos = 0.0*cm;  // y-direction (up-down of beam)
G4double NpolPolarimeter::HodoZPos = CHAnalyzerPos + 0.5*(HodoscopeZ + AnalyzerZ);

const G4Colour NpolPolarimeter::CopperColor  = G4Colour(0.0,1.0,0.0);
const G4Colour NpolPolarimeter::SteelColor  = G4Colour(0.5,1.0,0.5);
const G4Colour NpolPolarimeter::ScintColor  = G4Colour(0.0,0.5,0.25);

NpolPolarimeter::NpolPolarimeter() {
 
  copperAnalyzer = new NpolCopperAnalyzer();
  //gemDetectors = new NpolGEM();
  
}

NpolPolarimeter::~NpolPolarimeter() {}

G4String NpolPolarimeter::GetName() {
  return G4String("Polarimeter");
}

// There are "2" CH analyzer arrays.  This one is two 10cm by 10cm by 100cm long
// scintillator detectors in a vertical orientation of there long axis.
// They are placed on the outer edge of the acceptance with respect to the copper
// analyzer for the charge exchange process. 
void NpolPolarimeter::ConstructVertAnalyzer(G4LogicalVolume *motherLV) {

  G4double YPos = 0.0*m;  // y-direction (up-down of beam)
  G4double XPos = 0.0*m; // x-direction (left right of beam)
  G4double ZPos = CHAnalyzerPos;// CuAnalyzerPos + 22.925*cm; // // z-direction (in direction of beam)
  
  G4VSolid *VertAnalyzer = new G4Box("VertAnalyzer",vertAnalyzerX/2,vertAnalyzerY/2,vertAnalyzerZ/2);
  G4LogicalVolume *VertAnalyzerLV = new G4LogicalVolume(VertAnalyzer,
	  NpolMaterials::GetInstance()->GetMaterial("Scint"),"VertAnalyzerLV",0,0,0);
  
  G4AssemblyVolume *VertAnalyzerAV = new G4AssemblyVolume();
  G4RotationMatrix Ra, Rm; 
  G4ThreeVector Ta, Tm;
  G4Transform3D Tr;
  
  Ra.rotateX(0.0*deg);
  Ra.rotateY(0.0*deg);
  Ra.rotateZ(0.0*deg);
  // Place the first of two LV in the assembly
  Ta.setX(NpolCopperAnalyzer::CopperWidth/2-vertAnalyzerX/2); Ta.setY(0.0*cm); Ta.setZ(0.0*cm);
  Tr = G4Transform3D(Ra,Ta);
  VertAnalyzerAV->AddPlacedVolume(VertAnalyzerLV,Tr);
  // Place the second of two LV in the assembly
  Ta.setX(-(62.0*cm/2-vertAnalyzerX/2)); 
  Tr = G4Transform3D(Ra,Ta);
  VertAnalyzerAV->AddPlacedVolume(VertAnalyzerLV,Tr);
  // Place one imprint of the assembly in the world
  Tm.setX(XPos - ZPos*sin(NpolAng)); Tm.setY(YPos); Tm.setZ(ZPos*cos(NpolAng));
  Rm.rotateX(0.0*deg);
  Rm.rotateY(-NpolAng);
  Rm.rotateZ(0.0*deg);
  ImprintPlate(VertAnalyzerAV, motherLV, Tm, Rm);
  
  G4VisAttributes *TopVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  VertAnalyzerLV->SetVisAttributes(TopVisAtt);
}


//---------------------------
// Analyzer array based on University of Glasgow prototype highly segmented array.
// 4 x 8 scintillatory blocks each 4.0cm by 4.0cm by 25cm
//---------------------------
void NpolPolarimeter::ConstructGlasgowAnalyzer(G4LogicalVolume *motherLV) {
  
  G4double YPos = 0.0*m;  // y-direction (up-down of beam)
  G4double XPos = 0.0*m; // x-direction (left right of beam)
  G4double ZPos = CHAnalyzerPos; // z-direction (in direction of beam)
 
  G4VSolid *Analyzer = new G4Box("Analyzer",AnalyzerX/2,AnalyzerY/2,AnalyzerZ/2);
  G4LogicalVolume *AnalyzerLV = new G4LogicalVolume(Analyzer,
	  NpolMaterials::GetInstance()->GetMaterial("Scint"),"AnalyzerLV",0,0,0);
  
  G4AssemblyVolume *GlasgowAnalyzer = new G4AssemblyVolume();
  G4RotationMatrix Ra, Rm; 
  G4ThreeVector Ta, Tm;
  G4Transform3D Tr;

  Ra.rotateX(0.0*deg);
  Ra.rotateY(0.0*deg);
  Ra.rotateZ(0.0*deg);
  
  for(int i = 0; i < 8; i++) {
	for(int j = 0; j < 4; j++) {
	  Ta.setX(6.0*cm - 4.0*cm*j); Ta.setY(+14.0*cm-4.0*cm*i);
	  Ta.setZ(0.0*cm);
	  Tr = G4Transform3D(Ra,Ta);
	  GlasgowAnalyzer->AddPlacedVolume(AnalyzerLV,Tr);
	}
  }

  Tm.setX(XPos - ZPos*sin(NpolAng)); Tm.setY(YPos); Tm.setZ(ZPos*cos(NpolAng));
  Rm.rotateX(0.0*deg);
  Rm.rotateY(-NpolAng);
  Rm.rotateZ(0.0*deg);
  ImprintPlate(GlasgowAnalyzer, motherLV, Tm, Rm);
  
  G4VisAttributes *TopVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  AnalyzerLV->SetVisAttributes(TopVisAtt);
}

void NpolPolarimeter::ConstructHodoscopeArray(G4LogicalVolume *motherLV) {

  G4VSolid *LeftHodoscope = new G4Box("LeftHodoscope",LeftHodoscopeX/2,HodoscopeY/2,HodoscopeZ/2);
  G4VSolid *RightHodoscope = new G4Box("RightHodoscope",RightHodoscopeX/2,HodoscopeY/2,HodoscopeZ/2);
  G4LogicalVolume *LeftHodoscopeLV =
	new G4LogicalVolume(LeftHodoscope,NpolMaterials::GetInstance()->GetMaterial("Scint"), "LeftHodoscopeLV",0,0,0);
  G4LogicalVolume *RightHodoscopeLV =
	new G4LogicalVolume(RightHodoscope,NpolMaterials::GetInstance()->GetMaterial("Scint"), "RightHodoscopeLV",0,0,0);

  G4AssemblyVolume *LeftHodoscopeArray = new G4AssemblyVolume();
  G4AssemblyVolume *RightHodoscopeArray = new G4AssemblyVolume();
  G4RotationMatrix Ra, Rm; 
  G4ThreeVector Ta, Tm;
  G4Transform3D Tr;

  Ra.rotateX(0.0*deg);
  Ra.rotateY(0.0*deg);
  Ra.rotateZ(0.0*deg);
  
  for(int i = 0; i < 24; i++) {
	Ta.setX(0.0*cm); Ta.setY(+(11*HodoscopeY + 0.5*HodoscopeY) - HodoscopeY*i);
	Ta.setZ(0.0*cm);
	Tr = G4Transform3D(Ra,Ta);
	LeftHodoscopeArray->AddPlacedVolume(LeftHodoscopeLV,Tr);
  }

  Tm.setX(LeftHodoXPos*cos(NpolAng) + -HodoZPos*sin(NpolAng)); Tm.setY(HodoYPos);
  Tm.setZ(LeftHodoXPos*sin(NpolAng) + HodoZPos*cos(NpolAng));
  Rm.rotateX(0.0*deg);
  Rm.rotateY(-NpolAng);
  Rm.rotateZ(0.0*deg);
  ImprintPlate(LeftHodoscopeArray, motherLV, Tm, Rm); // left side
 
  for(int i = 0; i < 24; i++) {
	Ta.setX(0.0*cm); Ta.setY(+(11*HodoscopeY + 0.5*HodoscopeY) - HodoscopeY*i);
	Ta.setZ(0.0*cm);
	Tr = G4Transform3D(Ra,Ta);
	RightHodoscopeArray->AddPlacedVolume(RightHodoscopeLV,Tr);
  }
  
  Tm.setX(-RightHodoXPos*cos(NpolAng) + -HodoZPos*sin(NpolAng));
  Tm.setZ(-RightHodoXPos*sin(NpolAng) + HodoZPos*cos(NpolAng));
  ImprintPlate(RightHodoscopeArray, motherLV, Tm, Rm); // right side

  G4VisAttributes *HodoVisAtt= new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  LeftHodoscopeLV->SetVisAttributes(HodoVisAtt);
  RightHodoscopeLV->SetVisAttributes(HodoVisAtt);
}

void NpolPolarimeter::ConstructPolarimeterFluxTagger(G4LogicalVolume *motherLV){

  double width = 100*cm; double height = 150*cm; double thick = 0.1*cm;
  double xPos,zPos; double tagLocation = 380*cm;

  G4Box *NPOLTagger = new G4Box("NPOLTagger",width/2,height/2,thick/2);
  G4LogicalVolume *NPOLTaggerLV = new G4LogicalVolume(NPOLTagger,NpolMaterials::GetInstance()->GetMaterial("HardVacuum"),"NPOLTaggerLV",0,0,0);
  // G4VisAttributes *NPOLTaggerVisAtt = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2));
  //NPOLTaggerLV->SetVisAttributes(NPOLTaggerVisAtt);
  NPOLTaggerLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  xPos = -(tagLocation)*sin(NpolAng);
  zPos = +(tagLocation)*cos(NpolAng);

  PlaceRectangular(NPOLTaggerLV, motherLV, "NPOLTagger", xPos, 0.0*cm, zPos, 0*deg, -NpolAng, 0*deg);
}

void NpolPolarimeter::ConstructFakeGEM(G4LogicalVolume *motherLV){

  double infnGEMx = 40.0*cm, infnGEMy = 3*40.0*cm, infnGEMz = 0.250*cm;
  double uvaGEMx = 75.0*cm, uvaGEMy = 4*50.0*cm, uvaGEMz = 0.250*cm;

  double infnXpos = 0.0*cm ,infnYpos = 0.0*cm, infnZpos = 0.0*cm;
  double uvaXpos = 0.0*cm ,uvaYpos = 0.0*cm, uvaZpos = 0.0*cm;
  
  G4Box *infnGEMbox = new G4Box("infnGEMbox",infnGEMx/2,infnGEMy/2,infnGEMz/2);
  G4Box *uvaGEMbox = new G4Box("uvaGEMbox",uvaGEMx/2,uvaGEMy/2,uvaGEMz/2);

  G4LogicalVolume *infnGEMLV =
	new G4LogicalVolume(infnGEMbox,NpolMaterials::GetInstance()->GetMaterial("Scint")
						,"infnGEMLV",0,0,0);
  G4LogicalVolume *uvaGEMLV =
	new G4LogicalVolume(uvaGEMbox,NpolMaterials::GetInstance()->GetMaterial("Scint"),
						"uvaGEMLV",0,0,0);
  
  infnGEMLV->SetVisAttributes(ScintColor);
  uvaGEMLV->SetVisAttributes(ScintColor);
  
  G4RotationMatrix ra, rm;
  G4ThreeVector ta, tm;
  G4AssemblyVolume *infnModule = new G4AssemblyVolume();
  G4AssemblyVolume *uvaModule = new G4AssemblyVolume();
  G4AssemblyVolume *uvaModuleRear = new G4AssemblyVolume();
   
  ta.setX(0.); ta.setY(0.); ta.setZ(0.);
  infnModule->AddPlacedVolume(infnGEMLV,ta, &ra);
  uvaModule->AddPlacedVolume(uvaGEMLV,ta, &ra);
  ra.rotateY(90.0*deg);
  uvaModuleRear->AddPlacedVolume(uvaGEMLV, ta, &ra);
 
  rm.rotateX(0.0*deg);
  rm.rotateY(-NpolAng);
  rm.rotateZ(0.0*deg);

  // Two INFN GEMs in front of Cu analyzer
  infnZpos = CuAnalyzerPos + -39.53*cm;
  tm.setX(infnXpos-infnZpos*sin(NpolAng)); tm.setY(infnYpos); tm.setZ(infnZpos*cos(NpolAng));
  ImprintPlate(infnModule,motherLV, tm, rm);
  infnZpos = CuAnalyzerPos + -27.24*cm;
  tm.setX(infnXpos-infnZpos*sin(NpolAng)); tm.setY(infnYpos); tm.setZ(infnZpos*cos(NpolAng));
  ImprintPlate(infnModule,motherLV, tm, rm);

  // Two UVa GEMs in from of Cu analyzer
  uvaZpos = CuAnalyzerPos + -15.983*cm;
  tm.setX(-uvaZpos*sin(NpolAng)); tm.setY(uvaYpos); tm.setZ(uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModule,motherLV, tm, rm);
  uvaZpos = CuAnalyzerPos + -9.208*cm;
  tm.setX(-uvaZpos*sin(NpolAng)); tm.setY(uvaYpos); tm.setZ(uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModule,motherLV, tm, rm);

  // Four UVa GEMs behind Cu analyzer
  uvaZpos = CuAnalyzerPos + 9.207*cm;
  tm.setX(-uvaZpos*sin(NpolAng)); tm.setY(uvaYpos); tm.setZ(uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModule,motherLV, tm, rm);
  uvaZpos = CuAnalyzerPos + 16.547*cm;
  tm.setX(-uvaZpos*sin(NpolAng)); tm.setY(uvaYpos); tm.setZ(uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModule,motherLV, tm, rm);
  uvaZpos = CuAnalyzerPos + 29.303*cm;
  tm.setX(-uvaZpos*sin(NpolAng)); tm.setY(uvaYpos); tm.setZ(uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModule,motherLV, tm, rm);
  uvaZpos = CuAnalyzerPos + 42.073*cm;
  tm.setX(-uvaZpos*sin(NpolAng)); tm.setY(uvaYpos); tm.setZ(uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModule,motherLV, tm, rm);

  // Two UVa GEMs in from of hodoscopes on each side (4 total)
  uvaYpos = 0.0*cm; uvaZpos =  CHAnalyzerPos + 0.5*(HodoscopeZ + AnalyzerZ);
  uvaXpos = LeftHodoXPos - 9.63*cm; 
  tm.setX(uvaXpos*cos(NpolAng) + -uvaZpos*sin(NpolAng)); tm.setY(uvaYpos);
  tm.setZ(uvaXpos*sin(NpolAng) + uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModuleRear, motherLV, tm, rm);
 
  uvaXpos = LeftHodoXPos - 16.96*cm;
  tm.setX(uvaXpos*cos(NpolAng) + -uvaZpos*sin(NpolAng)); tm.setY(uvaYpos);
  tm.setZ(uvaXpos*sin(NpolAng) + uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModuleRear, motherLV, tm, rm);
  
  uvaXpos = -(RightHodoXPos - 10.98*cm);
  tm.setX(uvaXpos*cos(NpolAng) + -uvaZpos*sin(NpolAng)); tm.setY(uvaYpos);
  tm.setZ(uvaXpos*sin(NpolAng) + uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModuleRear, motherLV, tm, rm);
  
  uvaXpos = -(RightHodoXPos - 18.31*cm);
  tm.setX(uvaXpos*cos(NpolAng) + -uvaZpos*sin(NpolAng)); tm.setY(uvaYpos);
  tm.setZ(uvaXpos*sin(NpolAng) + uvaZpos*cos(NpolAng));
  ImprintPlate(uvaModuleRear, motherLV, tm, rm);
 
}



void NpolPolarimeter::Place(G4LogicalVolume *motherLV) {
  
  copperAnalyzer->Place(motherLV);
  ConstructVertAnalyzer(motherLV);
  //ConstructGlasgowAnalyzer(motherLV);
  ConstructFakeGEM(motherLV); // Just scintillator sheets in place of GEMs for tracking
  ConstructHodoscopeArray(motherLV);
  ConstructPolarimeterFluxTagger(motherLV);

  //gemDetectors->Place(motherLV); // "realistic" GEMs models
}


G4AssemblyVolume *NpolPolarimeter::MakePlate(G4LogicalVolume *detLV,	
	 G4int numDets, G4double TmX, G4double TmY, G4double TmZ,
	 G4double TmdX, G4double TmdY, G4double TmdZ) {
  
  G4AssemblyVolume *plate = new G4AssemblyVolume();
  
  // Translation and rotation of plate inside assembly
  G4RotationMatrix Ra; 
  G4ThreeVector Ta, Tm;
  G4Transform3D Tr;
  
  int i;
  for(i=0; i<numDets; i++) {
    Tm.setX(TmX-TmdX*i); Tm.setY(TmY-TmdY*i); Tm.setZ(TmZ-TmdZ*i);
    Tr = G4Transform3D(Ra,Tm);
    plate->AddPlacedVolume(detLV,Tr);
  }
  
  return plate;
}

void NpolPolarimeter::ImprintPlate(G4AssemblyVolume *plate, 
	G4LogicalVolume *motherLV, G4double TmX, G4double TmY, 
	G4double TmZ, G4double RmZ) {
  
  G4ThreeVector Tm;
  G4RotationMatrix Rm;
  G4Transform3D Tr;
  
  Tm.setX(TmX); Tm.setY(TmY); Tm.setZ(TmZ);
  Rm.rotateX(0.0*deg);
  Rm.rotateY(0.0*deg);
  Rm.rotateZ(RmZ);
  Tr = G4Transform3D(Rm,Tm);
  
  plate->MakeImprint(motherLV,Tr);
}

void NpolPolarimeter::ImprintPlate(G4AssemblyVolume *plate, G4LogicalVolume *motherLV, G4ThreeVector Tm, G4RotationMatrix Rm) {
  
  G4Transform3D Tr;
  Tr = G4Transform3D(Rm,Tm);
  
  plate->MakeImprint(motherLV,Tr);
}
