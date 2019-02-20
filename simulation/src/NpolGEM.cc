//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% NpolGEM.cc %%

// Constructor for GEM detectors;  lifted from TREK code provided by Michael Kohl (Hampton University).
// Created: William Tireman - November 2018

#include <stdio.h>
#include <cstdlib>
#include <math.h>

#include "G4PhysicalConstants.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4TwoVector.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4AssemblyVolume.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4String.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolGEM.hh"

const G4Colour NpolGEM::CopperColor  = G4Colour(0.0,1.0,0.0);
const G4Colour NpolGEM::MylarColor   = G4Colour(0.0,1.0,1.0);
const G4Colour NpolGEM::NpolGasColor = G4Colour(1.0,0.0,1.0);
const G4Colour NpolGEM::KaptonColor  = G4Colour(0.0,1.0,1.0);

NpolGEM::NpolGEM(){

  ConstructINFNGEM();
  ConstructUVaGEM();
}

NpolGEM::~NpolGEM(){}

G4String NpolGEM::GetName() {
  return G4String("GEM Detector");
}


void NpolGEM::ConstructINFNGEM(){
  
  //GEM Assembly Volume    
  assemblyINFNGEM = new G4AssemblyVolume();
  BuildModule(40*cm,50*cm); // in centimeters
  GEMLV();
  ConstructDetector(assemblyINFNGEM,3, 40.0*cm, 40.0*cm); // in centimeters
  
}

void NpolGEM::ConstructUVaGEM(){
  assemblyUVaGEM = new G4AssemblyVolume();
  BuildModule(50*cm,60*cm); // in centimeters
  GEMLV();
  ConstructDetector(assemblyUVaGEM,4, 75.0*cm, 50.0*cm);
  
}

void NpolGEM::BuildModule(G4double height, G4double width){

  GEMFoil_x = 0.0025*mm;
  GEMFoil_y = height; 
  GEMFoil_z = width; 
  
  GEMKapton_x = 0.025*mm;
  GEMKapton_y = height; 
  GEMKapton_z = width; 
  
  GEMMylar_x = 0.010*mm;
  
  GEMgapReadout_x = 3.00*mm;
  
  GEMgap_x = 1.5*mm;
  GEMgap_y = height; 
  GEMgap_z = width; 
  
  //GEM geometry definition
  GEMFoil = 
    new G4Box("GEMFoil",GEMFoil_x/2, GEMFoil_y/2, GEMFoil_z/2);
  
  GEMKapton = 
    new G4Box("GEMFoilKapton", GEMKapton_x/2, GEMKapton_y/2, GEMKapton_z/2);
  
  GEMgap = 
    new G4Box("GEMgap", GEMgap_x/2, GEMgap_y/2, GEMgap_z/2);
  
  GEMMylar = 
    new G4Box("GEMFoil", GEMMylar_x/2, GEMKapton_y/2, GEMKapton_z/2);
  
  GEMdriftCathode = 
    new G4Box("GEMFoil", GEMgap_x/2, GEMgap_y/2, GEMgap_z/2);
  
  GEMReadout = 
    new G4Box("GEMreadout", GEMKapton_x/2, GEMKapton_y/2, GEMKapton_z/2);
  
  GEMdriftR = 
    new G4Box("GEMFoil", GEMgapReadout_x/2, GEMKapton_y/2, GEMKapton_z/2);
  
}

void NpolGEM::ConstructDetector(G4AssemblyVolume *assemblyGEM,G4double numModules,G4double vertPosStart,G4double vertShift) {
  
  //Rotation and translation for GEMAssVol
  G4RotationMatrix ra;
  G4ThreeVector ta;
  
  //Translation and Rotation of GEM Assembly in World Volume
  G4ThreeVector tm;
  G4RotationMatrix rm;
  G4cout<<"<----------Rotation of GEMAssVol------------->\n";
  
  G4double gemPos = 0.00*mm;
  G4double y;
  //GEM foils
  for(G4int i=0; i < numModules; i ++){
	phi = 0.*deg;//*(i);
	y = vertPosStart - (vertShift * i);
	rm.rotateZ(phi);
	
	//Fill Assembly with GEM components
	ta.setX(gemPos); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMMylarLV(), ta, &ra);
	
	ta.setX(gemPos + 1.525*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMdriftCathLV(), ta, &ra);
	
	ta.setX(gemPos + 3.030*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMFoilLV(), ta, &ra);
	
	ta.setX(gemPos + 3.055*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMKaptonLV(), ta, &ra);
	
	ta.setX(gemPos + 3.085*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMFoilLV(), ta, &ra); 
	
	ta.setX(gemPos + 4.625*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMgapLV(), ta, &ra); 
	
	ta.setX(gemPos + 6.5*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMFoilLV(), ta, &ra);
	
	ta.setX(gemPos + 6.55*mm); ta.setY(y); ta.setZ(0.); 
	assemblyGEM->AddPlacedVolume(GEMKaptonLV(), ta, &ra);
	
 	ta.setX(gemPos + 6.7*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMFoilLV(), ta, &ra);
	
	ta.setX(gemPos + 9.*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMgapLV(), ta, &ra);
	
	ta.setX(gemPos + 11.*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMFoilLV(), ta, &ra);
	
	ta.setX(gemPos + 11.125*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMKaptonLV(), ta, &ra);
	
	ta.setX(gemPos + 11.25*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMFoilLV(), ta, &ra);
	
	ta.setX(gemPos + 14.25*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMdriftReadLV(), ta, &ra);
	
	ta.setX(gemPos + 17.25*mm); ta.setY(y); ta.setZ(0.);
	assemblyGEM->AddPlacedVolume(GEMReadoutLV(), ta, &ra);
	
  }
}


void NpolGEM::GEMLV() {
  //GEM Logical Volume construction
  GEMFoil_log = new G4LogicalVolume(GEMFoil,NpolMaterials::GetInstance()->GetMaterial("Cu"),GEMFoil_name);
  
  GEMKapton_log = new G4LogicalVolume(GEMKapton,NpolMaterials::GetInstance()->GetMaterial("Kapton"),GEMKapton_name); 
  
  GEMMylar_log = new G4LogicalVolume(GEMMylar,NpolMaterials::GetInstance()->GetMaterial("Mylar"),GEMFoil_name);  
  
  GEMdriftCathode_log = new G4LogicalVolume(GEMdriftCathode,NpolMaterials::GetInstance()->GetMaterial("NpolGas"),GEMFoil_name);  
  
  GEMReadout_log = new G4LogicalVolume(GEMReadout,NpolMaterials::GetInstance()->GetMaterial("Cu"),GEMKapton_name);
  GEMdriftReadout_log = new G4LogicalVolume(GEMdriftR,NpolMaterials::GetInstance()->GetMaterial("NpolGas"),GEMKapton_name);
  GEMgap_log =  new G4LogicalVolume(GEMgap,NpolMaterials::GetInstance()->GetMaterial("NpolGas"),GEMgap_name);
  
  //Visualization Attributes
  GEMMylar_log         -> SetVisAttributes(MylarColor);
  GEMFoil_log          -> SetVisAttributes(CopperColor);
  GEMKapton_log        -> SetVisAttributes(KaptonColor);
  GEMgap_log           -> SetVisAttributes(NpolGasColor);
  GEMdriftCathode_log  -> SetVisAttributes(NpolGasColor);
  GEMReadout_log       -> SetVisAttributes(CopperColor);
  GEMdriftReadout_log  -> SetVisAttributes(NpolGasColor);
}


void NpolGEM::Place(G4LogicalVolume *motherLV){
  
  //Translation and Rotation of GEM Assembly in World Volume
  G4ThreeVector tm;  G4RotationMatrix RotINFN, RotUVaFront, RotUVaRear;

  // Place the GEM detectors before the Cu analyzer
  RotINFN.rotateX(0.0*deg); RotINFN.rotateY(90.0*deg); RotINFN.rotateZ(0.0*deg);
  tm.setX(0.0*cm); tm.setY(0.0*cm); tm.setZ(-116.0*cm);
  assemblyINFNGEM->MakeImprint(motherLV, tm, &RotINFN, 0, false);
  tm.setZ(-108.0*cm);
  assemblyINFNGEM->MakeImprint(motherLV, tm, &RotINFN, 0, false);
  
  RotUVaFront.rotateX(0.0*deg); RotUVaFront.rotateY(90.0*deg); RotUVaFront.rotateZ(0.0*deg);
  tm.setZ(-93.0*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaFront, 0, false);
  tm.setZ(-100.5*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaFront, 0, false);
 

  // Place the GEM detectors after the Cu analyzer
  tm.setZ(-70.5*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaFront, 0, false);
  tm.setZ(-63.0*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaFront, 0, false);
  tm.setZ(-54.0*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaFront, 0, false);
  tm.setZ(-46.0*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaFront, 0, false);

  // Place the large angle GEM detectors for (n,p) elastic scattering
  RotUVaRear.rotateX(0.0*deg); RotUVaRear.rotateY(0.0*deg); RotUVaRear.rotateZ(0.0*deg);
  tm.setX(-90.0*cm); tm.setZ(10.0*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaRear, 0, false);
  tm.setX(-97.5*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaRear, 0, false);
  RotUVaRear.rotateY(+180.0*deg);
  tm.setX(+90.0*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaRear, 0, false);
  tm.setX(+97.5*cm);
  assemblyUVaGEM->MakeImprint(motherLV, tm, &RotUVaRear, 0, false);
}
