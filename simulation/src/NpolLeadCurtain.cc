//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	         *
//******************************************************************

// %% NpolLeadCurtain.cc %%

// Lead Curtain construction file
// Created: William Tireman - May 2017

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolLeadCurtain.hh"
#include "NpolDipole2.hh"
#include "NpolParticleFluxTagger.hh"
#include "NpolPolarimeter.hh"

G4double NpolLeadCurtain::leadThickness = 15.0*cm;  // thickness of the lead curtain
G4double NpolLeadCurtain::PosLead =  NpolDipole2::PosD2 - NpolDipole2::ClampOffSetF - NpolDipole2::fieldClampThick/2 - leadThickness/2 - 0.5*cm;

NpolLeadCurtain::NpolLeadCurtain() {
   ConstructLeadCurtain();
}

NpolLeadCurtain::~NpolLeadCurtain() {}

G4String NpolLeadCurtain::GetName() {
  return G4String("Lead Curtain");
}

// Construct a lead shield for in front of the polarimeter
void NpolLeadCurtain::ConstructLeadCurtain(){

  G4double xlen = 2*(PosLead)*tan(NpolParticleFluxTagger::horAngle/2);// + 10.0*cm;//80*cm; 
  G4double ylen = 2*(PosLead)*tan(NpolParticleFluxTagger::vertAngle/2);// + 10.0*cm; //80*cm;

  G4Box *LeadCurtain = new G4Box("LeadCurtain",xlen/2,ylen/2,leadThickness/2);
  LeadCurtainLV = new G4LogicalVolume(LeadCurtain,NpolMaterials::GetInstance()->GetMaterial("Pb"),"LeadCurtainLV",0,0,0);
  G4VisAttributes *LeadCurtainVisAtt = new G4VisAttributes(G4Colour(1.0, 0.9, 0.6));
  LeadCurtainLV->SetVisAttributes(LeadCurtainVisAtt);
}


void NpolLeadCurtain::Place(G4LogicalVolume *motherLV) {
 
  PlaceCylindrical(LeadCurtainLV, motherLV, "LeadCurtain", PosLead,-NpolPolarimeter::NpolAng, 0);

}
