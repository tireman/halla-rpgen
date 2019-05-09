//********************************************************************
//* License and Disclaimer: From GEANT Collaboration                 *
//*                                                                  *
//* The  Geant4 software  is  copyright of the Copyright Holders  of *
//* the Geant4 Collaboration.  It is provided  under  the terms  and *
//* conditions of the Geant4 Software License,  included in the file *
//* LICENSE and available at  http://cern.ch/geant4/license .  These *
//* include a list of copyright holders.     		      	*
//********************************************************************

// %% NpolSBSBeamline.cc %%

// SBS Beamline Constructor file.  
// Created: William Tireman - November 2018
// Heavily modified from SBS simulation (g4sbs)

#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Cons.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolSBSBeamline.hh"
#include "NpolPolarimeter.hh"
#include "NpolDipole2.hh"

G4int NpolSBSBeamline::fBeamlineConf = 3; //GMn
G4String NpolSBSBeamline::TargType = "kLD2"; // Opts:"kH2", "k3He", "kNeutTarg"
G4String NpolSBSBeamline::ExpType = "kGMn";  // Opts: "kGEp" "kGMn"
G4int NpolSBSBeamline::ScattChamberFlag = 1;
G4int NpolSBSBeamline::LeadOption = 1;
G4double NpolSBSBeamline::f48D48depth = 1219.2*mm;
G4double NpolSBSBeamline::f48D48width = 2324.1*mm;
G4double NpolSBSBeamline::f48D48height = 3721.1*mm;
G4double NpolSBSBeamline::f48D48ang  = NpolPolarimeter::NpolAng; //24.7*deg; //39.4*deg;
G4double NpolSBSBeamline::f48D48dist = NpolDipole2::PosD2; //2.8*m;
G4double NpolSBSBeamline::f48D48_fieldclamp_config = 2;// 0 = No field clamps. 2 = GEp (default). 1 = BigBite experiments:
G4double NpolSBSBeamline::beamheight = 10.0*12*2.54*cm; // 10 feet off the ground
G4double NpolSBSBeamline::gapheight= 700.*mm;
G4double NpolSBSBeamline::ent_len = 10*m;
G4double NpolSBSBeamline::swallrad_inner = 1.041/2.0*m;
G4double NpolSBSBeamline::swallrad = 1.143*m/2;
G4double NpolSBSBeamline::ent_rin = 31.75*mm;
G4double NpolSBSBeamline::ent_rou = ent_rin+0.120*mm;
G4double NpolSBSBeamline::shieldblock3_height = (f48D48depth - gapheight)/2;
G4double NpolSBSBeamline::leadstart = 290*cm;
G4double NpolSBSBeamline::leadend   = 435*cm;
G4double NpolSBSBeamline::magleadlen = leadend-leadstart;

NpolMaterials *Mat = NpolMaterials::GetInstance();

NpolSBSBeamline::NpolSBSBeamline() {
  
  //ConstructGMnBeamline();
  ConstructGEpLead();
  //ConstructGMnLead();
  
}

NpolSBSBeamline::~NpolSBSBeamline() {}

G4String NpolSBSBeamline::GetName() {
  return G4String("SBS Beamline");
}



void NpolSBSBeamline::ConstructGEpLead(){
  
  G4VisAttributes *lead_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  
  G4double inch = 2.54*cm;
  G4double TargetCenter_zoffset = 6.50*inch;
  
  G4double z_outer_magnetic = 182.33*cm - TargetCenter_zoffset;
  
  G4double zstart_lead1 = 170.0*cm;
  G4double z_formed_bellows = 133.2*cm - TargetCenter_zoffset;
  G4double zstop_lead1 = z_formed_bellows + 75.0*inch;
  
  //G4cout << "zmag, zstart, zstop = " << z_outer_magnetic << ", " << zstart_lead1 << ", " << zstop_lead1 << G4endl;
  
  G4double Rin1 = 9.3*cm;
  G4double Rout1 = Rin1 + 5.0*cm;
  G4double Rin2 = Rin1 + (zstop_lead1 - zstart_lead1)*tan(1.5*deg );
  G4double Rout2 = Rin2 + 5.0*cm;
  
  G4Cons *leadcone1 = new G4Cons("leadcone1", Rin1, Rout1, Rin2, Rout2, (zstop_lead1-zstart_lead1)/2.0, 90.0*deg, 180.0*deg );
  
  G4double width = 0.5*f48D48width;
  G4double angle = f48D48ang;
  G4double dist = f48D48dist;
  G4double depth = f48D48depth;
  
  //We want to subtract the overlap between the cone and the SBS magnet.
  G4Box *cutbox1 = new G4Box( "cutbox1", width/2.0, Rout2+cm, depth/2.0 );
  G4Box *slottemp = new G4Box( "slottemp", width/2.0, 15.5*cm, depth/2.0 + cm );
  G4RotationMatrix *rot_temp = new G4RotationMatrix;
  rot_temp->rotateY(angle);
  
  G4SubtractionSolid *boxwithslot = new G4SubtractionSolid( "boxwithslot", cutbox1, slottemp, 0, G4ThreeVector(35.0*cm,0,0) );
  //G4LogicalVolume *boxwithslot_log = new G4LogicalVolume( boxwithslot, Mat->GetMaterial("Air"), "boxwithslot_log" );
  
  G4double Rbox = dist + depth/2.0;
  
  G4ThreeVector pos_box_withslot( -Rbox*sin(angle) + width/2.0*cos(angle), 0, Rbox*cos(angle) + width/2.0*sin(angle) );
  
  //new G4PVPlacement( rot_temp, pos_box_withslot, boxwithslot_log, "boxwithslot_phys", worldlog, false, 0 );
  
  G4ThreeVector pos(0,0,0.5*(zstart_lead1+zstop_lead1));
  
  G4ThreeVector posrel_boxwithslot = pos_box_withslot - pos;
  
  G4SubtractionSolid *leadcone1_cut = new G4SubtractionSolid("leadcone1_cut", leadcone1, boxwithslot, rot_temp, posrel_boxwithslot );
  
  leadcone1_log = new G4LogicalVolume(leadcone1_cut, Mat->GetMaterial("Pb"), "leadcone1_log" );
  
  leadcone1_log->SetVisAttributes( lead_visatt );
  
  
  //new G4PVPlacement( 0, pos, leadcone1_log, "leadcone1_phys", worldlog, false, 0 );
  
  G4double zsections[3] = {z_outer_magnetic + 74.0*inch,
						   201.632*inch - TargetCenter_zoffset,
						   207.144*inch - TargetCenter_zoffset + 40.0*inch };
  G4double Rin_sections[3] = { 15.0*inch/2.0 + (zsections[0]-zsections[1])*tan(1.5*deg),
							   15.0*inch/2.0,
							   15.0*inch/2.0 };
  G4double Rout_sections[3] = {Rin_sections[0] + 5.*cm,
							   Rin_sections[1] + 5.*cm,
							   Rin_sections[2] + 5.*cm };
  
  G4Polycone *leadshield2 = new G4Polycone( "leadshield2", 90.0*deg, 180.0*deg, 3, zsections, Rin_sections, Rout_sections );
  leadshield2_log = new G4LogicalVolume( leadshield2, Mat->GetMaterial("Pb"), "leadshield2_log" );
  
  leadshield2_log->SetVisAttributes( lead_visatt );
  
  //new G4PVPlacement( 0, G4ThreeVector(), leadshield2_log, "leadshield2_phys", worldlog, false, 0 );
  
  ////// New geometry with vertical wall(s): 
  
  G4ThreeVector zaxis_temp( -sin(16.9*deg), 0.0, cos(16.9*deg) );
  G4ThreeVector yaxis_temp( 0,1,0);
  G4ThreeVector xaxis_temp = (yaxis_temp.cross(zaxis_temp)).unit();
  
  G4ThreeVector frontcorner_pos = 1.6*m*zaxis_temp;
  
  G4double zstart_lead_wall1 = z_outer_magnetic + 15*cm;
  //G4double zstop_lead_wall1 = zstart_lead_wall1 + 1.25*m;
  
  G4Box *lead_wall1 = new G4Box("lead_wall1", 5.0*cm/2, 30.5*cm/2, 1.25*m/2 ); 
  lead_wall1_log = new G4LogicalVolume( lead_wall1, Mat->GetMaterial("Pb"), "lead_wall1_log" );
  
  G4double xtemp = -( 5.5*inch/2.0 + 1.5*inch + (1.25/2.0+0.15)*m*tan(1.5*deg) + 2.5*cm/cos(1.5*deg) );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( 1.5*deg );
  
  G4cout << "Lead wall A (x,y,z) = (" << xtemp/cm << ", " << 0.0 << ", " << (zstart_lead_wall1 + 0.5*1.25*m)/cm << ")" << G4endl;
  
  lead_wall1_log->SetVisAttributes( lead_visatt );
  
  G4double zstart_lead_wall2 = z_formed_bellows + 76.09*inch + 1.71*inch + 15.75*inch + 1.0*inch;
  G4double zstop_lead_wall2 = 207.144*inch - TargetCenter_zoffset + 38.8*inch;  /// W.T. changed, original 40*inch
  
  G4cout << "Lead wall B zstart - zstop = " << (zstop_lead_wall2 - zstart_lead_wall2)/cm << G4endl;
  
  G4double zpos_lead_wall2 = 0.5*(zstart_lead_wall2 + zstop_lead_wall2 );
  //we want x position to have x = 
  G4double xpos_lead_wall2 = -(8.0*inch + 2.5*cm + (zpos_lead_wall2 - 201.632*inch + TargetCenter_zoffset )*tan(1.5*deg));
  
  G4Box *lead_wall2 = new G4Box("lead_wall2", 5.08*cm/2, 42.00*inch/2, 2.4384*m/2 /*0.5*(zstop_lead_wall2 - zstart_lead_wall2)*/ ); // 18-April-2019 <-- W.T. reduced sizes to those in CAD drawing from Feb. 2019

  lead_wall2_log = new G4LogicalVolume( lead_wall2, Mat->GetMaterial("Pb"), "lead_wall2_log" );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( 1.5*deg );
  
  G4cout << "Lead wall B (x,y,z) = (" << xpos_lead_wall2/cm << ", " << 0.0 << ", " << zpos_lead_wall2/cm << ")" << G4endl;
  
  lead_wall2_log->SetVisAttributes( lead_visatt );
  
}


void NpolSBSBeamline::ConstructCommonExitBeamline(G4LogicalVolume *worldlog) {
 //Define visualization attributes here:
  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.1, 0.5, 0.9 ) );
  //Vacuum_visatt->SetForceSolid(true);
  Vacuum_visatt->SetVisibility(false);
  G4VisAttributes *SteelColor = new G4VisAttributes( G4Colour( 0.75, 0.75, 0.75 ) );
  G4VisAttributes *CopperColor = new G4VisAttributes( G4Colour( 0.7, 0.3, 0.3 ) );
  
  G4double inch = 2.54*cm;
  
  G4double TargetCenter_zoffset = 6.50*inch;
  
  G4double z_formed_bellows = 52.440*inch - TargetCenter_zoffset; //relative to "target center"? or "origin"?
  G4double z_spool_piece = 58.44*inch - TargetCenter_zoffset;
  if(fBeamlineConf>2)z_spool_piece = 27.903*inch;
  G4double z_conic_vacline_weldment = 62.8*inch - TargetCenter_zoffset;
  G4double z_outer_magnetic = 71.782*inch - TargetCenter_zoffset;
  G4double z_inner_magnetic = 73.782*inch - TargetCenter_zoffset;
  G4double z_welded_bellows = 201.632*inch - TargetCenter_zoffset;
  
  G4double X=0.0, Y=0.0, Z=0.0;
  G4ThreeVector zero(0.0, 0.0, 0.0);
  
  // Conic vacuum line weldment upstream flange:

  G4double Rin, Rout, Thick;
  G4double Rin1, Rout1, Rin2, Rout2;
  Rin = 3.517*inch/2.0;
  Rout = 6.75*inch/2.0;
  
  Thick = 0.84*inch;
  Rin1 = Rin - Thick/2.0*tan( 1.5*deg );
  Rout1 = Rout;
  Rin2 = Rin + Thick/2.0*tan( 1.5*deg );
  Rout2 = Rout;
  
  //G4Cons *CVLW_Flange1 = new G4Cons( "CVLW_Flange1", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange1 = new G4Tubs("CVLW_Flange1", Rin, Rout, Thick/2.0, 0.0, twopi );
  
  // Fill with vacuum
  Rin = 0.0;
  Rout = 3.517*inch/2.0;
  
  Rout1 = Rout - Thick/2.0 * tan( 1.5*deg );
  Rout2 = Rout + Thick/2.0 * tan( 1.5*deg );
  
  //G4Cons *CVLW_Flange1_vac = new G4Cons("CVLW_Flange1_vac", Rin, Rout1, Rin, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange1_vac = new G4Tubs("CVLW_Flange1_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
  
  G4LogicalVolume *CVLW_Flange1_log = new G4LogicalVolume( CVLW_Flange1, Mat->GetMaterial("SSteel"), "CVLW_Flange1_log" );
  G4LogicalVolume *CVLW_Flange1_vac_log = new G4LogicalVolume( CVLW_Flange1_vac, Mat->GetMaterial("Vacuum"), "CVLW_Flange1_vac_log" );
  
  CVLW_Flange1_log->SetVisAttributes( SteelColor );
  CVLW_Flange1_vac_log->SetVisAttributes( Vacuum_visatt );
  
  // Then place the vacuum inside the Iron Tube
  Z = z_conic_vacline_weldment + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_vac_log, "CVLW_Flange1_vac_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_log, "CVLW_Flange1_phys", worldlog, false, 0 );
  
  

  //conic vacuum line weldment:
  //Thick = 3.50*m;
  Thick = 138.83*inch - 1.12*inch - 0.84*inch;
  Rin1 = 3.517*inch/2.0;
  Rout1 = Rin1 + 0.125*inch;
  Rin2 = 10.734*inch/2.0;
  Rout2 = Rin2 + 0.125*inch;
  
  G4Cons *CVLW = new G4Cons( "CVLW", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );   
  // Fill with vacuum				
  
  Rin1 = 0.0;
  Rout1 = 3.517*inch/2.0;
  Rin2 = 0.0;
  Rout2 = 10.734*inch/2.0;
  
  G4Cons *CVLW_vac = new G4Cons( "CVLW_vac", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  // Convert into logical volumes
  G4LogicalVolume *CVLW_log = new G4LogicalVolume( CVLW, Mat->GetMaterial("SSteel"), "CVLW_log" );
  G4LogicalVolume *CVLW_vac_log = new G4LogicalVolume( CVLW_vac, Mat->GetMaterial("Vacuum"), "CVLW_vac_log" );
  
  CVLW_log->SetVisAttributes( SteelColor );
  CVLW_vac_log->SetVisAttributes( Vacuum_visatt );
  // Then place the vacuum inside the Iron Cone
  //Z = (159.51 + 2.13)*cm + pDz;
  Z = z_conic_vacline_weldment + 0.84*inch + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_vac_log, "CVLW_vac_phys", worldlog, false, 0);
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_log, "CVLW_phys", worldlog, false, 0 );
  
  // Flange 2:
  Rin = 10.734/2.0*inch;
  Rout = 14.0/2.0*inch;
  Thick = 1.12*inch;
  Rin1 = Rin - Thick/2.0 * tan(1.5*deg);
  Rout1 = Rout;
  Rin2 = Rin + Thick/2.0 * tan(1.5*deg);
  Rout2 = Rout;
  //G4Cons *CVLW_Flange2 = new G4Cons("CVLW_Flange2", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange2 = new G4Tubs("CVLW_Flange2", Rin, Rout, Thick/2.0, 0.0, twopi );
  // Fill with vacuum
  Rin = 0.0;
  Rout1 = Rin1;
  Rout2 = Rin2;
  Rout = 10.734*inch/2.0;
  
  //G4Cons *CVLW_Flange2_vac = new G4Cons("CVLW_Flange2_vac", Rin, Rout1, Rin, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange2_vac = new G4Tubs( "CVLW_Flange2_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
  // Convert into logical volumes
  // G4LogicalVolume *FLN2_log = new G4LogicalVolume( FLN2_tube, Mat->GetMaterial("Fe"), "FLN2_log" );
  // G4LogicalVolume *FVL2_log = new G4LogicalVolume( FVL2_tube, Mat->GetMaterial("Vacuum"), "FVL2_log");
  G4LogicalVolume *CVLW_Flange2_log = new G4LogicalVolume( CVLW_Flange2, Mat->GetMaterial("SSteel"), "CVLW_Flange2_log" );
  G4LogicalVolume *CVLW_Flange2_vac_log = new G4LogicalVolume( CVLW_Flange2_vac, Mat->GetMaterial("Vacuum"), "CVLW_Flange2_vac_log" );
  
  CVLW_Flange2_log->SetVisAttributes(SteelColor );
  CVLW_Flange2_vac_log->SetVisAttributes( Vacuum_visatt );

  // Then place the vacuum inside the Iron Tube
  Z = z_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_vac_log, "CVLW_Flange2_vac_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_log, "CVLW_Flange2_phys", worldlog, false, 0 );

  //Next: "Welded bellows"
  G4double dz_welded_bellows = 207.144*inch - z_welded_bellows - TargetCenter_zoffset; // = =5.512 inches
  
  Rin = 11.750/2.0*inch;
  Rout = 14.0/2.0*inch;
  Thick = 1.12*inch;
  
  G4Tubs *WB_Flange = new G4Tubs( "WB_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Flange_log = new G4LogicalVolume( WB_Flange, Mat->GetMaterial("SSteel"), "WB_Flange_log" );
  
  WB_Flange_log->SetVisAttributes( SteelColor );

  Z = z_welded_bellows + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange1_phys", worldlog, false, 0 );
  
  Z = z_welded_bellows + dz_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange2_phys", worldlog, false, 1 );
  
  Rout = Rin + 0.125*inch;
  Thick = dz_welded_bellows - 2*1.12*inch;
  G4Tubs *WB_Bellows = new G4Tubs( "WB_Bellows", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Bellows_log = new G4LogicalVolume(WB_Bellows, Mat->GetMaterial("SSteel"), "WB_Bellows_log" );
  
  WB_Bellows_log->SetVisAttributes( SteelColor );
  
  Z = z_welded_bellows + 1.12*inch + Thick/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Bellows_log, "WB_Bellows_phys", worldlog, false, 0 );
  
  Rin = 0.0;
  Rout = 11.750/2.0*inch;
  Thick = dz_welded_bellows;
  G4Tubs *WB_Vacuum = new G4Tubs( "WB_Vacuum", Rin, Rout, Thick/2.0, 0.0, twopi );
  
  G4LogicalVolume *WB_Vacuum_log = new G4LogicalVolume(WB_Vacuum, Mat->GetMaterial("Vacuum"), "WB_Vacuum_log" );
  
  WB_Vacuum_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_welded_bellows + dz_welded_bellows/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Vacuum_log, "WB_Vacuum_phys", worldlog, false, 0 );

  // Here a bellow and we assign wall of 0.03 cm
  // tRmin = (0.5*27.62)*cm;
  // tRmax = (0.5*27.62 + 0.03)*cm;
  // tDzz  = 0.5*(4.237*2.54 - 2.84)*cm;
  // G4Tubs *TBL8_tube = new G4Tubs("TBL8_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // // Fill with vacuum
  // tRmin = 0.0*cm;
  // tRmax = (0.5*27.62)*cm;
  // G4Tubs *TVL8_tube = new G4Tubs("TVL8_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // // Convert into logical volumes
  // G4LogicalVolume *TBL8_log = new G4LogicalVolume( TBL8_tube, Mat->GetMaterial("Fe"), "TBL8_log" );
  // G4LogicalVolume *TVL8_log = new G4LogicalVolume( TVL8_tube, Mat->GetMaterial("Vacuum"), "TVL8_log" );
  // // Then place the vacuum inside the Iron Tube
  // Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54*0.5)*cm;
  // new G4PVPlacement( 0, zero, TVL8_log, "Bellow_Vac", TBL8_log, false, 0 );
  // new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL8_log, "Bellow_Iron", worldlog, false, 0 );

  //return;
  
  // EXTEND VACUUM LINE by using Maduka geometry
  // ============================================
  G4double tRmin, tRmax, tDzz, pDz, pRmax1, pRmax2, tSPhi, tDphi, pRmin1, pRmin2;
  tSPhi = 0.0;
  tDphi = twopi;
  
  tRmin = 0.5*12.0*2.54*cm; 
  tRmax = 13.0*2.54*0.5*cm;
  tDzz  = 0.5*41.0*2.54*cm;
  G4Tubs *TBL9_tube = new G4Tubs("TBL9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Fill with vacuum
  tRmin = 0.0*cm;
  tRmax = 0.5*12.0*2.54*cm;
  G4Tubs *TVL9_tube = new G4Tubs("TVL9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Convert into logical volumes
  G4LogicalVolume *TBL9_log = new G4LogicalVolume( TBL9_tube, Mat->GetMaterial("Al"), "TBL9_log" );
  G4LogicalVolume *TVL9_log = new G4LogicalVolume( TVL9_tube, Mat->GetMaterial("Vacuum"), "TVL9_log" );

  TBL9_log->SetVisAttributes( AlColor );
  TVL9_log->SetVisAttributes( Vacuum_visatt );
  
  // Then place the vacuum inside the Al Tube
  //Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54*0.5)*cm;
  Z = 207.144*inch + tDzz - TargetCenter_zoffset;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TVL9_log, "Extended_Vac1", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL9_log, "Extended_Al1", worldlog, false, 0 );

  tRmin = 0.5*24.0*2.54*cm;
  tRmax = 25.0*2.54*0.5*cm;
  tDzz  = 0.5*217.0*2.54*cm;
  G4Tubs *TML9_tube = new G4Tubs( "TML9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Fill with vacuum
  tRmin = 0.0*cm;
  tRmax = 0.5*24.0*2.54*cm;
  G4Tubs *TMV9_tube = new G4Tubs("TMV9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Convert into logical volumes
  G4LogicalVolume *TML9_log = new G4LogicalVolume( TML9_tube, Mat->GetMaterial("Al"), "TML9_log" );
  G4LogicalVolume *TMV9_log = new G4LogicalVolume( TMV9_tube, Mat->GetMaterial("Vacuum"), "TMV9_log" );

  TML9_log->SetVisAttributes( AlColor );
  TMV9_log->SetVisAttributes( Vacuum_visatt );
  // Then place vacuum inside of Al tube
  //Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54 + 0.5*217.0*2.54)*cm;
  Z = 207.144*inch + 41.0*inch + tDzz - TargetCenter_zoffset;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TMV9_log, "Extended_Vac2", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TML9_log, "Extended_Al2", worldlog, false, 0 );

  // For CPU speed, extend vacuum all the way to the edge of the "world" volume, so that we don't track beam electrons in air beyond interesting region.
  G4double Zstop = 50.0*m;
  G4double Zstart = Z + tDzz;
  G4double Zwidth = (Zstop-Zstart);
  G4Tubs *FakeVacuumExtension = new G4Tubs( "FakeVacuumExtension", tRmin, tRmax, Zwidth/2.0, tSPhi, tDphi );
  G4LogicalVolume *FakeVacuumExtension_log = new G4LogicalVolume( FakeVacuumExtension, Mat->GetMaterial("Vacuum"), "FakeVacuumExtension_log" );
  FakeVacuumExtension_log->SetVisAttributes( Vacuum_visatt );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0.5*(Zstop+Zstart)), FakeVacuumExtension_log, "FakeVacuumExtension_phys", worldlog,false,0);

  //-----------------------------------------------------
  //       magnetic tubes
    
  // Inner and outer Magnetic shieldings: 
  // arrays of variables to parameterize the different geometries with the beamline config flag
  G4int Ndivs = 1;// number of segments with shielding
  //G4double Rin_array[6];// radii for inner shielding elements
  //G4double Zin_array[6];// z for inner shielding elements
  //G4int Nrings_out[3];// number of outer elements per segments
  //G4double Rout_array[6];// radii for outer shielding elements
  //G4double Zout_array[6];// z for outer shielding elements
  std::vector<G4double> Rin_array;// radii for inner shielding elements
  std::vector<G4double> Zin_array;// z for inner shielding elements
  std::vector<G4int> Nrings_out;// number of outer elements per segments
  std::vector<G4double> Rout_array;// radii for inner shielding elements
  std::vector<G4double> Zout_array;// z for inner shielding elements

  G4double OMthick = 1.625*inch;
  G4double OMspace = 0.375*inch;

  switch(fBeamlineConf){
  case(1):// reminder: beamline config 1 = GEp
    Ndivs = 2;
    
    Rin_array.push_back( 4.354*inch/2.0 );
    Rin_array.push_back( 6.848*inch/2.0 );
    Rin_array.push_back( 8.230*inch/2.0 );
    Rin_array.push_back( 10.971*inch/2.0 );
    
    Zin_array.push_back( z_inner_magnetic );
    Zin_array.push_back( z_inner_magnetic + 47.625*inch );
    Zin_array.push_back( z_inner_magnetic + 74.00*inch );
    Zin_array.push_back( z_inner_magnetic + (74.00 + 52.347)*inch );
    
    Nrings_out.push_back( 26 );
    Nrings_out.push_back( 27 );
    
    Rout_array.push_back( 5.5*inch/2.0 );
    Rout_array.push_back( 8.178*inch/2.0 );
    Rout_array.push_back( 9.349*inch/2.0 );
    Rout_array.push_back( 12.156*inch/2.0 );
    
    Zout_array.push_back( z_outer_magnetic );
    Zout_array.push_back( z_outer_magnetic + 26.0*OMthick + 25.0*OMspace );
    Zout_array.push_back( z_outer_magnetic + 74.0*inch );
    Zout_array.push_back( z_outer_magnetic + 74.0*inch + 27.0*OMthick + 26.0*OMspace );
    
    //G4double Rin_array_tmp1[6] = {4.354*inch/2.0, 6.848*inch/2.0, 8.230*inch/2.0, 10.971*inch/2.0, 0.0, 0.0};
    //G4double Zin_array_tmp1[6] = {z_inner_magnetic, z_inner_magnetic + 47.625*inch, z_inner_magnetic + 74.00*inch, z_inner_magnetic + (74.00 + 52.347)*inch, 0.0, 0.0};
    //G4int Nrings_out_tmp1[3] = {26, 27, 0};
    //G4double Rout_array_tmp1[6] = {5.5*inch/2.0, 8.178*inch/2.0, 9.349*inch/2.0, 12.156*inch/2.0, 0.0, 0.0};
    //G4double Zout_array_tmp1[6] = {z_outer_magnetic, z_outer_magnetic + 26.0*OMthick + 25.0*OMspace, z_outer_magnetic + 74.0*inch, z_outer_magnetic + 74.0*inch + 27.0*OMthick + 26.0*OMspace, 0.0, 0.0};
    break;
  // case(2):// reminder: beamline config 2 = GEn, SIDIS
  // Ndivs = 3;
  //   break;
  case(3):// reminder: beamline config 3 = GMn, all Q^2
    Ndivs = 3;
    
    Rin_array.push_back( 3.774*inch/2.0 );
    Rin_array.push_back( 4.307*inch/2.0 );
    Rin_array.push_back( 5.264*inch/2.0 );
    Rin_array.push_back( 7.905*inch/2.0 );
    Rin_array.push_back( 9.310*inch/2.0 );
    Rin_array.push_back( 10.930*inch/2.0 );
    
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 +11.62 - 1.625)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 +11.62 + 14.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 - 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 138.83 - 1.12 - 1.14)*inch );
    
    Nrings_out.push_back( 6 );
    Nrings_out.push_back( 27 );
    Nrings_out.push_back( 17 );
    
    Rout_array.push_back( (3.774/2+0.38)*inch );
    Rout_array.push_back( (4.392/2+0.38)*inch );
    Rout_array.push_back( (5.158/2+0.38)*inch );
    Rout_array.push_back( (8.012/2+0.38)*inch );
    Rout_array.push_back( (9.203/2+0.38)*inch );
    Rout_array.push_back( (10.923/2+0.38)*inch );
    
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );
    
    //G4double Rin_array_tmp3[6] = {3.774*inch/2.0, 4.307*inch/2.0, 5.264*inch/2.0, 7.905*inch/2.0, 9.310*inch/2.0, 10.930*inch/2.0};
    //G4double Zin_array_tmp3[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62 - 1.625)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 - 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 138.83 - 1.12 - 1.14)*inch};
    //G4int Nrings_out_tmp3[3] = {6, 27, 17};
    //G4double Rout_array_tmp3[6] = {(3.774/2+0.38)*inch, (4.392/2+0.38)*inch, (5.158/2+0.38)*inch, (8.012/2+0.38)*inch, (9.203/2+0.38)*inch, (10.923/2+0.38)*inch};
    //G4double Zout_array_tmp3[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38)*inch, z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch};
    break;
  case(4):// reminder: beamline config 4 = GMn, Q^2 = 13.5 GeV^2 (+calibrations)
    Ndivs = 2;
    
    Rin_array.push_back( 3.774*inch/2.0 );
    Rin_array.push_back( 6.096*inch/2.0 );
    Rin_array.push_back( 7.074*inch/2.0 );
    Rin_array.push_back( 9.609*inch/2.0 );
    
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 - 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62 - 2.0)*inch );
    
    Nrings_out.push_back( 23 );
    Nrings_out.push_back( 26 );
    
    Rout_array.push_back( (3.774/2+0.38)*inch );
    Rout_array.push_back( (6.202/2+0.38)*inch );
    Rout_array.push_back( (6.968/2+0.38)*inch );
    Rout_array.push_back( (9.715/2+0.38)*inch );
    
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch  );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62)*inch );
    
    //G4double Rin_array_tmp4[6] = {3.774*inch/2.0, 6.096*inch/2.0, 7.074*inch/2.0, 9.609*inch/2.0, 0.0, 0.0};
    //G4double Zin_array_tmp4[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 - 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62 - 2.0)*inch, 0.0, 0.0};
    //G4int Nrings_out_tmp4[3] = {23, 26, 0};
    //G4double Rout_array_tmp4[6] = {(3.774/2+0.38)*inch, (6.202/2+0.38)*inch, (6.968/2+0.38)*inch, (9.715/2+0.38)*inch, 0.0, 0.0};
    //G4double Zout_array_tmp4[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62)*inch, 0.0, 0.0};    
    break;
  default:// default: all elements built
    Ndivs = 1;
    Rin_array.push_back( 3.774*inch/2.0 );
    Rin_array.push_back( 10.923*inch/2.0 );

    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );

    Nrings_out.push_back( 68 );
    
    Rout_array.push_back( (3.774/2+0.38)*inch );
    Rout_array.push_back( (10.923/2+0.38)*inch );

    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );
    
    //G4double Rin_array_tmp0[6] = {3.774*inch/2.0, 10.923*inch/2.0, 0.0, 0.0, 0.0, 0.0};
    //G4double Zin_array_tmp0[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch, 0.0, 0.0, 0.0, 0.0};
    //G4int Nrings_out_tmp0[3] = {68, 0, 0};
    //G4double Rout_array_tmp0[6] = {(3.774/2+0.38)*inch, (10.923/2+0.38)*inch, 0.0, 0.0, 0.0, 0.0};
    //G4double Zout_array_tmp0[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch, 0.0, 0.0, 0.0, 0.0};
    break;
  }
  
  // Building beamline shielding:
  for(G4int i = 0; i<Ndivs; i++){
    // Building beamline shielding: inner elements
    Rin1 = Rin_array[2*i];
    Rout1 = Rin1 + 0.25*inch;
    Rin2 = Rin_array[2*i+1];
    Rout2 = Rin2 + 0.25*inch;
    Thick = Zin_array[2*i+1]-Zin_array[2*i];
	
    char cname[100];
    sprintf(cname,"IM_%d", i);
    G4String name = cname;
     
    G4Cons *IM_ = new G4Cons( name, Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
    name += "_log";
    G4LogicalVolume *IM__log = new G4LogicalVolume( IM_, Mat->GetMaterial("Fe"), name );
    
    IM__log->SetVisAttributes( ironColor );
        
    Z = (Zin_array[2*i]+Zin_array[2*i+1])/2.0;
    
    name = cname;
    name += "_phys";
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), IM__log, name, worldlog, false, 0 );
    
    // Building beamline shielding: outer elements
    G4double zmin = Zout_array[2*i];
    G4double zmax = Zout_array[2*i+1];
    
    G4double Rin_min = Rout_array[2*i];
    G4double Rin_max = Rout_array[2*i+1];
    for( G4int j=0; j<Nrings_out[i]; j++ ){
      char cname2[100];
      sprintf(cname2,"OM_%d_ring%d", i, j);
      G4String name2 = cname2;
      
      G4double zstart = zmin + j*(OMthick + OMspace);
      G4double zstop = zstart + OMthick;
      G4double Rin_start = Rin_min + (zstart-zmin)/(zmax-zmin)*(Rin_max - Rin_min);
      G4double Rout_start = Rin_start + 0.5*inch;
      G4double Rin_stop = Rin_min + (zstop-zmin)/(zmax-zmin)*(Rin_max - Rin_min);
      G4double Rout_stop = Rin_stop + 0.5*inch;
   
      G4Cons *ring = new G4Cons( name2,Rin_start, Rout_start, Rin_stop, Rout_stop, OMthick/2.0, 0.0, twopi );
	
      name2 += "_log";
      G4LogicalVolume *ring_log = new G4LogicalVolume( ring, Mat->GetMaterial("Fe"), name2 );
      
      ring_log->SetVisAttributes( ironColor );
      
      name2 = cname2;
      name2 += "_phys";
      
      Z = 0.5*(zstart + zstop);
      new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name2, worldlog, false, 0 );
    }
  }
 
  if(fBeamlineConf != 2){
    G4double dz_spool_piece = z_conic_vacline_weldment - z_spool_piece;
    
    //Make Spool piece vacuum:
    Rout = 3.76*inch/2.0;
    Rin = 0.0;
    Thick = dz_spool_piece;
    
    G4Tubs *SpoolPiece_vac = new G4Tubs("SpoolPiece_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *SpoolPiece_vac_log = new G4LogicalVolume( SpoolPiece_vac, Mat->GetMaterial("Vacuum"), "SpoolPiece_vac_log" );
  
    SpoolPiece_vac_log->SetVisAttributes( Vacuum_visatt );
    
    Z = z_spool_piece + Thick/2.0;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), SpoolPiece_vac_log, "SpoolPiece_vac_phys", worldlog, false, 0 );
    
    Rin = 3.76*inch/2.0;
    Rout = 6.00*inch/2.0;
    Thick = 0.84*inch;
    
    G4Tubs *SpoolPiece_Flange1 = new G4Tubs("SpoolPiece_Flange1", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *SpoolPiece_Flange1_log = new G4LogicalVolume( SpoolPiece_Flange1, Mat->GetMaterial("SSteel"), "SpoolPiece_Flange1_log" );
  
    SpoolPiece_Flange1_log->SetVisAttributes( SteelColor );
    
    Z = z_spool_piece + Thick/2.0;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), SpoolPiece_Flange1_log, "SpoolPiece_Flange1_phys", worldlog, false, 0 );
    
    Rout = 6.75*inch/2.0;
    Thick = 0.84*inch;
    
    G4Tubs *SpoolPiece_Flange2 = new G4Tubs("SpoolPiece_Flange2", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *SpoolPiece_Flange2_log = new G4LogicalVolume( SpoolPiece_Flange2, Mat->GetMaterial("SSteel"), "SpoolPiece_Flange2_log" );
    
    SpoolPiece_Flange2_log->SetVisAttributes( SteelColor );
    
    Z = z_conic_vacline_weldment - Thick/2.0;
    
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), SpoolPiece_Flange2_log, "SpoolPiece_Flange2_phys", worldlog, false, 0 );
    
    Rout = 4.0*inch/2.0;
    Thick = dz_spool_piece - 2.0*0.84*inch;
    
    G4Tubs *SpoolPiece_tube = new G4Tubs("SpoolPiece_tube", Rin, Rout, Thick/2.0, 0.0, twopi );
    
    G4LogicalVolume *SpoolPiece_tube_log = new G4LogicalVolume( SpoolPiece_tube, Mat->GetMaterial("SSteel"), "SpoolPiece_tube_log" );
    
    SpoolPiece_tube_log->SetVisAttributes( SteelColor );
    
    Z = z_spool_piece + dz_spool_piece/2.0;
    
    new G4PVPlacement( 0,  G4ThreeVector( X, Y, Z ), SpoolPiece_tube_log, "SpoolPiece_tube_phys", worldlog, false, 0 );
  
	if(fBeamlineConf > 1){
	  Rin1 = 4.00*inch/2.0+0.02*inch;
	  Rout1 = Rin1 + 0.25*inch;
	  //Thick = dz_spool_piece - 2.0*0.84*inch;
	  
	  G4Tubs *IM0 = new G4Tubs( "IM0", Rin1, Rout1, Thick/2.0, 0.0, twopi );
	  G4LogicalVolume *IM0_log = new G4LogicalVolume( IM0, Mat->GetMaterial("Fe"), "IM0_log" );
	  
	  IM0_log->SetVisAttributes( ironColor );
	  
	  //Z = z_spool_piece + dz_spool_piece/2.0;
	  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), IM0_log, "IM0_phys", worldlog, false, 0 );
	
	  G4double zmin = z_spool_piece+0.84*inch+0.44*inch;
	  G4double zmax = zmin + 13.0*OMthick + 12.0*OMspace;
	
	  G4double Rin_min = 5.3*inch/2.0;
	  for( G4int i=0; i<13; i++ ){
		char cname3[100];
		sprintf(cname3,"OM0_ring%d", i);
		G4String name3 = cname3;
	  
		G4double zstart = zmin + i*(OMthick + OMspace);
		G4double zstop = zstart + OMthick;
	  
		G4Tubs *ring = new G4Tubs( name3,Rin_min, Rin_min+0.5*inch, OMthick/2.0, 0.0, twopi );
	  
		name3 += "_log";
		G4LogicalVolume *ring_log = new G4LogicalVolume( ring, Mat->GetMaterial("Fe"), name3 );
	 
		ring_log->SetVisAttributes( ironColor );
   
		name3 = cname3;
		name3 += "_phys";
	  
		Z = 0.5*(zstart + zstop);
		new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name3, worldlog, false, 0 );
	  }
	}
  }

  if(fBeamlineConf == 1){
    //Last but not least: formed bellows! defer to tomorrow...
   
    G4double dz_formed_bellows = 6.00*inch;
    Rin = 0.0;
    Rout = 3.81*inch/2.0;
    Thick = dz_formed_bellows;
    //define vacuum volume for formed bellows
    G4Tubs *FormedBellows_vac = new G4Tubs("FormedBellows_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *FormedBellows_vac_log = new G4LogicalVolume( FormedBellows_vac, Mat->GetMaterial("Vacuum"), "FormedBellows_vac_log" );
  
    FormedBellows_vac_log->SetVisAttributes( Vacuum_visatt );
    
    Z = z_formed_bellows + Thick/2.0;
    
    new G4PVPlacement(  0,  G4ThreeVector( X, Y, Z ), FormedBellows_vac_log, "FormedBellows_vac_phys", worldlog, false, 0 );
  
    Rin = 3.81*inch/2.0;
    Rout = 6.00*inch/2.0;
    Thick = 0.84*inch;
  
    //Flanges for formed bellows:
    G4Tubs *FormedBellows_Flange = new G4Tubs("FormedBellows_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *FormedBellows_Flange_log = new G4LogicalVolume( FormedBellows_Flange, Mat->GetMaterial("SSteel"), "FormedBellows_Flange_log" );
	
    FormedBellows_Flange_log->SetVisAttributes( SteelColor );
    
    Z = z_formed_bellows + Thick/2.0;
    
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_Flange_log, "FormedBellows_Flange1_phys", worldlog, false, 0 );
    
    Z = z_formed_bellows + dz_formed_bellows - Thick/2.0;
    
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_Flange_log, "FormedBellows_Flange2_phys", worldlog, false, 1 );
    
    //Tube for formed bellows:
    
    Rout = Rin + 0.125*inch; //This is just a guess!!
    Thick = dz_formed_bellows - 2.0*0.84*inch;
    
    G4Tubs *FormedBellows_tube = new G4Tubs( "FormedBellows_tube", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *FormedBellows_tube_log = new G4LogicalVolume( FormedBellows_tube, Mat->GetMaterial("SSteel"), "FormedBellows_tube_log" );
	
    FormedBellows_tube_log->SetVisAttributes( SteelColor );
	
    Z = z_formed_bellows + dz_formed_bellows/2.0;
	
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_tube_log, "FormedBellows_tube_phys", worldlog, false, 0 );
	
    //Two more "Iron" tubes to connect Snout to "formed bellows"
    G4double dz_iron_tubes = z_formed_bellows - 49.56*inch + TargetCenter_zoffset;
	
    Thick = dz_iron_tubes/2.0;
    Rin = 5.0*cm;
    Rout = 7.0*cm;
	
    G4Tubs *IronTube1 = new G4Tubs("IronTube1", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *IronTube1_log = new G4LogicalVolume( IronTube1, Mat->GetMaterial("Fe"), "IronTube1_log" );
    IronTube1_log->SetVisAttributes( ironColor );
  
    G4Tubs *IronTube1_vac = new G4Tubs("IronTube1_vac", 0.0, Rin, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *IronTube1_vac_log = new G4LogicalVolume( IronTube1_vac, Mat->GetMaterial("Vacuum"), "IronTube1_vac_log" );
  
    IronTube1_vac_log->SetVisAttributes( Vacuum_visatt );

    Z = 49.56*inch + Thick/2.0 - TargetCenter_zoffset;

    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube1_log, "IronTube1_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube1_vac_log, "IronTube1_vac_phys", worldlog, false, 0 );

    Rin = 2.415*inch;
    Rout = 2.5*inch;
  
    G4Tubs *IronTube2 = new G4Tubs("IronTube2", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4Tubs *IronTube2_vac = new G4Tubs("IronTube2_vac", 0.0, Rin, Thick/2.0, 0.0, twopi );

    G4LogicalVolume *IronTube2_log = new G4LogicalVolume( IronTube2, Mat->GetMaterial("Fe"), "IronTube2_log" );
    G4LogicalVolume *IronTube2_vac_log = new G4LogicalVolume( IronTube2_vac, Mat->GetMaterial("Vacuum"), "IronTube2_vac_log" );
  
    IronTube2_log->SetVisAttributes(ironColor);
    IronTube2_vac_log->SetVisAttributes(Vacuum_visatt);

    Z += Thick;
  
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube2_log, "IronTube2_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube2_vac_log, "IronTube2_vac_phys", worldlog, false, 0 );
  }

  //Next, corrector magnets:
  // Define some dimensions that are going to be useful to define the distances
  G4double UpstreamCoilThickY = 1.68*inch;
  G4double UpstreamCoilThickX = 3.46*inch;
  //G4double UpstreamCoilWidth = 3.46*inch;
  G4double UpstreamCoilHeight = 8.17*inch;
  G4double UpstreamCoilDepth = 6.60*inch;
  G4double UpstreamCoilWidth = 7.56*inch;
 
  G4double YokeTopPiece_Width = 15.04*inch;
  G4double YokeTopPiece_Height = 3.94*inch;
  G4double YokeTopPiece_Depth = 6.30*inch;
  
  G4double YokeLeftPiece_Width = 2.76*inch;
  G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;
  
  G4double YokeRightNotchAngle = 18.43*deg;
  G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );
  
  G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  G4double DownstreamYokeDepth = 15.75*inch;

  G4double DS_coil_depth = 8.91*inch;
  G4double DS_coil_height = 12.04*inch;
  G4double DS_coil_ThickX = 2.90*inch;
  G4double DS_coil_ThickY = 1.68*inch;
  
  // Z Array to change easily z values with beamline configuration;
  // Right now, it looks like X and Y do NOT need to change depending on the configuration; only Z does
  std::vector<G4double> z_Magnets_array;
  switch(fBeamlineConf){
  case(1):// reminder: beamline config 1 = GEp
    z_Magnets_array.push_back( z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
    z_Magnets_array.push_back( z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //z_Magnets_array = {z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY, z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0, z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0, z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0};
     break;
  // case(2):// reminder: beamline config 2 = GEn, SIDIS ?
  //   Z = z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  //   break;
  case(3):// reminder: beamline config 3 = GMn 
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 15.94)*inch + UpstreamCoilDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 15.94 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 85.78)*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 85.78 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //z_Magnets_array = {z_conic_vacline_weldment + (0.84 + 0.14 + 15.94)*inch + UpstreamCoilDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 15.94 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 85.78)*inch + DownstreamYokeDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 85.78 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0};    
    break;
  case(4):// reminder: beamline config 3 = GMn, Q^2 = 13.5 GeV^2
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 49.47)*inch + UpstreamCoilDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 49.47 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0  );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 118.34)*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 118.34 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
  
    //z_Magnets_array = {z_conic_vacline_weldment + (0.84 + 0.14 + 49.47)*inch + UpstreamCoilDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 49.47 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 118.34)*inch + DownstreamYokeDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 118.34 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0};
    break;
  default:
    //sent back: who cares...
    z_Magnets_array.push_back( -10.0*m + 0.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
    z_Magnets_array.push_back( -10.0*m + 2.3*inch + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( -10.0*m + 15.24*inch + 1.71*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( -10.0*m + 15.24*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    //z_Magnets_array = {z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY, z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0, z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0, z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0};
    break;
  }

  G4Box *UpstreamCoil_outer = new G4Box("UpstreamCoil_outer", UpstreamCoilThickX/2.0, (UpstreamCoilHeight+2.0*UpstreamCoilThickY)/2.0, (UpstreamCoilDepth + 2.0*UpstreamCoilThickY)/2.0 );
  G4Box *UpstreamCoil_inner = new G4Box("UpstreamCoil_inner", UpstreamCoilThickX/2.0 + cm, UpstreamCoilHeight/2.0, UpstreamCoilDepth/2.0 );

  G4SubtractionSolid *UpstreamCoil = new G4SubtractionSolid( "UpstreamCoil", UpstreamCoil_outer, UpstreamCoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *UpstreamCoil_log = new G4LogicalVolume(UpstreamCoil, Mat->GetMaterial("Cu"), "UpstreamCoil_log" );

  UpstreamCoil_log->SetVisAttributes( CopperColor );

  Z = z_Magnets_array[0];//z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  X = (UpstreamCoilWidth+UpstreamCoilThickX)/2.0;
  Y = 0.0;

  //two placements of upstream coil:
  
  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_right", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_left", worldlog, false, 1 );

  G4double UpstreamPoleDepth = 6.3*inch;
  G4double UpstreamPoleWidth = 4.02*inch;
  G4double UpstreamPoleHeight = 7.87*inch;
  //Next, make poles:
  G4Box *UpstreamPole = new G4Box( "UpstreamPole", UpstreamPoleWidth/2.0, UpstreamPoleHeight/2.0, UpstreamPoleDepth/2.0 );
  G4LogicalVolume *UpstreamPole_log = new G4LogicalVolume( UpstreamPole, Mat->GetMaterial("Fe"), "UpstreamPole_log" );
  UpstreamPole_log->SetVisAttributes( ironColor );
  //two placements of upstream poles:

  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_right", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_left", worldlog, false, 1 );

  //Next, make surrounding yoke:
  // G4double YokeTopPiece_Width = 15.04*inch;
  // G4double YokeTopPiece_Height = 3.94*inch;
  // G4double YokeTopPiece_Depth = 6.30*inch;

  G4Box *YokeTopPiece = new G4Box("YokeTopPiece", YokeTopPiece_Width/2.0, YokeTopPiece_Height/2.0, YokeTopPiece_Depth/2.0 );
  G4LogicalVolume *YokeTopPiece_log = new G4LogicalVolume( YokeTopPiece, Mat->GetMaterial("Fe"), "YokeTopPiece_log" );

  YokeTopPiece_log->SetVisAttributes( ironColor );
  
  X = 0.0;
  Y = (11.81*inch + YokeTopPiece_Height)/2.0;

  //two placements of yoke top piece (top and bottom symmetric):
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeTopPiece_log, "UpstreamYokeTop_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X,-Y,Z), YokeTopPiece_log, "UpstreamYokeBottom_phys", worldlog, false, 1 );

  // G4double YokeLeftPiece_Width = 2.76*inch;
  // G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  // G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;
  
  G4Box *YokeLeftPiece = new G4Box("YokeLeftPiece", YokeLeftPiece_Width/2.0, YokeLeftPiece_Height/2.0, YokeLeftPiece_Depth/2.0 );
  G4LogicalVolume *YokeLeftPiece_log = new G4LogicalVolume( YokeLeftPiece, Mat->GetMaterial("Fe"), "YokeLeftPiece_log" );
  YokeLeftPiece_log->SetVisAttributes(ironColor );

  X = 7.52*inch + YokeLeftPiece_Width/2.0;
  Y = 0.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeLeftPiece_log, "UpstreamYokeLeftPiece_phys", worldlog, false, 0 );

  // G4double YokeRightNotchAngle = 18.43*deg;
  // G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  // G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  // G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );

  //I think this is correct:
  G4Trap *YokeRight_trap = new G4Trap( "YokeRight_trap", YokeRightZFinal/2.0, atan( (YokeRightWidthFinal-YokeRightWidthInitial)/2.0/YokeRightZFinal ), 180.0*deg,
				       YokeLeftPiece_Height/2.0, YokeRightWidthInitial/2.0, YokeRightWidthInitial/2.0, 0.0,
				       YokeLeftPiece_Height/2.0, YokeRightWidthFinal/2.0, YokeRightWidthFinal/2.0, 0.0 ); 

  G4Box *YokeRight_box = new G4Box( "YokeRight_box", YokeRightWidthFinal/2.0, YokeLeftPiece_Height/2.0, 0.39*inch/2.0 );

  X = 0.5*(YokeRightWidthFinal/2.0 + YokeRightWidthInitial/2.0) - YokeRightWidthFinal/2.0;
  
  G4UnionSolid *YokeRightPiece = new G4UnionSolid("YokeRightPiece", YokeRight_trap, YokeRight_box, 0, G4ThreeVector( X, 0, (YokeRightZFinal+0.39*inch)/2.0 ) );
  G4LogicalVolume *YokeRightPiece_log = new G4LogicalVolume(YokeRightPiece, Mat->GetMaterial("Fe"), "YokeRightPiece_log" );

  YokeRightPiece_log->SetVisAttributes(ironColor);

  X = -7.52*inch - 0.5*(YokeRightWidthFinal/2.0 + YokeRightWidthInitial/2.0);
  Y = 0.0;
  Z = z_Magnets_array[1];//z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeRightPiece_log, "UpstreamYokeRightPiece_phys", worldlog, false, 0 );

  //Downstream Corrector:
  // G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  // G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  // G4double DownstreamYokeDepth = 15.75*inch;

  G4double DownstreamYokeGapWidth = 17.58*inch;
  G4double DownstreamYokeGapHeight = 20.16*inch;
  G4Box *DownstreamYoke_box = new G4Box("DownstreamYoke_box", DownstreamTotalWidth/2.0, DownstreamTotalHeight/2.0, DownstreamYokeDepth/2.0 );
  G4Box *DownstreamYoke_gap = new G4Box("DownstreamYoke_gap", DownstreamYokeGapWidth/2.0, DownstreamYokeGapHeight/2.0, DownstreamYokeDepth/2.0+cm );
  G4SubtractionSolid *DownstreamYoke = new G4SubtractionSolid( "DownstreamYoke", DownstreamYoke_box, DownstreamYoke_gap, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DownstreamYoke_log = new G4LogicalVolume( DownstreamYoke, Mat->GetMaterial("Fe"), "DownstreamYoke_log" );

  DownstreamYoke_log->SetVisAttributes( ironColor );

  X = 0.0; Y = 0.0;
  Z = z_Magnets_array[2];//z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DownstreamYoke_log, "DownstreamYoke_phys", worldlog, false, 0 );

  // G4double DS_coil_depth = 8.91*inch;
  // G4double DS_coil_height = 12.04*inch;
  // G4double DS_coil_ThickX = 2.90*inch;
  // G4double DS_coil_ThickY = 1.68*inch;
  
  G4Box *DS_coil_outer = new G4Box( "DS_coil_outer", DS_coil_ThickX/2.0, (DS_coil_height + 2.0*DS_coil_ThickY)/2.0, (DS_coil_depth + 2.0*DS_coil_ThickY)/2.0 );
  G4Box *DS_coil_inner = new G4Box( "DS_coil_inner", DS_coil_ThickX/2.0+cm, DS_coil_height/2.0, DS_coil_depth/2.0 );

  G4SubtractionSolid *DS_coil = new G4SubtractionSolid( "DS_coil", DS_coil_outer, DS_coil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DS_coil_log = new G4LogicalVolume( DS_coil, Mat->GetMaterial("Cu"), "DS_coil_log" );
  DS_coil_log->SetVisAttributes(CopperColor );
  
  X = 11.67*inch/2.0 + DS_coil_ThickX/2.0;
  Y = 0.0;
  Z = z_Magnets_array[3];//z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DS_coil_log, "DS_coil_phys_left", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z), DS_coil_log, "DS_coil_phys_right", worldlog, false, 1 );

  //Now just need poles:
  G4double DSpole_depth = 8.76*inch;
  G4double DSpole_width = (17.58-11.00)*inch/2.0;
  G4double DSpole_height = 11.81*inch;

  G4Box *DSpole = new G4Box("DSpole", DSpole_width/2.0, DSpole_height/2.0, DSpole_depth/2.0 );
  G4LogicalVolume *DSpole_log = new G4LogicalVolume(DSpole, Mat->GetMaterial("Fe"), "DSpole_log" );

  DSpole_log->SetVisAttributes(ironColor);
  
  X = (17.58+11.00)*inch/4.0;
  Y = 0.0;
  //two placements of poles:
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DSpole_log, "DSpole_phys_left", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z), DSpole_log, "DSpole_phys_right", worldlog, false, 1 );
  
}

  /*if(fDetCon->fBLneutronDet){//TO-DO: set the possibility to deactivate it.
    // EFuchey: 2018/05/29: add a small dummy detector to study the neutron production by the shielding.
    // 
    // double x_blndet[8] = {0.3*m, 0.6*m, 1.6*m, 0.0*m, -0.7*m, 1.0*m, 2.2*m, -2.2*m};
    // double y_blndet[8] = {0.0*m, 0.0*m, 0.0*m, 0.6*m,  0.0*m, 0.0*m, 0.0*m,  0.0*m};
    // double z_blndet[8] = {1.5*m, 2.5*m, 5.0*m, 2.5*m, +0.7*m, 0.0*m, 2.2*m,  2.2*m};
    //
    // G4double ElecX = 2.0*cm;
    // G4double ElecY = 2.0*cm;
    // G4double ElecZ = 2.0*cm;

    double x_blndet = 3.0*m;
    double y_blndet = 0.0*m;
    double z_blndet = 2.5*m;
    
    G4double ElecX = 5.0*cm;
    G4double ElecY = 100.0*cm;
    G4double ElecZ = 100.0*cm;
    
    G4Box *Electronics = new G4Box( "Electronics" , ElecX/2.0, ElecY/2.0, ElecZ/2.0);
    G4LogicalVolume *Electronics_log = new G4LogicalVolume( Electronics , Mat->GetMaterial("Silicon"), "Electronics_log" );
    
    G4String GEMElectronicsname = "BLneutronDet";
    G4String  GEMElectronicscollname = "BLneutronDet";
    G4SBSCalSD *GEMElecSD = NULL;
    
    switch(fDetCon->fExpType){
    case(kGEp):
      GEMElectronicsname += "GEp";
      GEMElectronicscollname += "GEp";
      break;
    case(kNeutronExp):// GMn
      GEMElectronicsname += "GMn";
      GEMElectronicscollname += "GMn";
      break;
    case(kGEnRP):// GEnRP
      GEMElectronicsname += "GMn";
      GEMElectronicscollname += "GMn";
      break;
    default:
      
      break;
    }
    
    //for(int i_blndet = 0; i_blndet<8; i_blndet++){
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(GEMElectronicsname) )){
      G4cout << "Adding GEM electronics Sensitive Detector to SDman..." << G4endl;
      GEMElecSD = new G4SBSCalSD( GEMElectronicsname, GEMElectronicscollname );
      fDetCon->fSDman->AddNewDetector(GEMElecSD);
      (fDetCon->SDlist).insert(GEMElectronicsname);
      fDetCon->SDtype[GEMElectronicsname] = kCAL;
      (GEMElecSD->detmap).depth = 0;
    }
    Electronics_log->SetSensitiveDetector( GEMElecSD );
    
    if( (fDetCon->StepLimiterList).find( GEMElectronicsname ) != (fDetCon->StepLimiterList).end() ){
      Electronics_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }
    
    // Place the electronics in our hut:
    // new G4PVPlacement( 0, G4ThreeVector(0.0, -ShieldMotherY/2.0 + GPlateY2 + ElecY/2.0, ShieldMotherZ/2.0 - GPlateZ1 - ElecZ/2.0),
    // 		     Electronics_log, "Electronics", ShieldLog, false, 0);
    // new G4PVPlacement( 0, G4ThreeVector(x_blndet[i_blndet], y_blndet[i_blndet], z_blndet[i_blndet]),
    // 		       Electronics_log, "GMn_Electronics", worldlog, false, i_blndet);
    new G4PVPlacement( 0, G4ThreeVector(x_blndet, y_blndet, z_blndet),
		       Electronics_log, "GMn_Electronics", worldlog, false, 0);
  }*/
  //}

  // VISUALS
  
  // CVLW_Flange1_log->SetVisAttributes( ironColor );
  // CVLW_log->SetVisAttributes( ironColor );
  // CVLW_Flange2_log->SetVisAttributes( ironColor );
  // WB_Flange_log->SetVisAttributes( ironColor );
  // WB_Bellows_log->SetVisAttributes( ironColor );
  // //TBL8_log->SetVisAttributes( ironColor );
  // TBM1_log->SetVisAttributes( ironColor );
  // TBM2_log->SetVisAttributes( ironColor );
  // TBM3_log->SetVisAttributes( ironColor );
  // TBM4_log->SetVisAttributes( ironColor );
  // TBT1_log->SetVisAttributes( ironColor );
  // TBT2_log->SetVisAttributes( ironColor );

  
  // TBL9_log->SetVisAttributes( AlColor );
  // TML9_log->SetVisAttributes( AlColor );

  // // Vacuum
  // FVL1_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL2_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL3_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL5_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL6_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL7_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TVB1_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TVL8_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TVL9_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TMV9_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TTV1_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TTV2_log->SetVisAttributes( G4VisAttributes::Invisible );





// This is the beam line for GMn
void NpolSBSBeamline::ConstructGMnBeamline(G4LogicalVolume *worldlog){
  double swallrad = 1.143*m/2;
  double swallrad_inner = 1.041/2.0*m; 
  //EFuchey: 2017/02/14: change parameters for Standard scat chamber:
  double sc_entbeampipeflange_dist = 25.375*2.54*cm;// entrance pipe flange distance from hall center
  double sc_exbeampipeflange_dist = 27.903*2.54*cm;// exit pipe flange distance from hall center
  
  // Stainless
  G4double ent_len = 10*m;
  //ent_len = ent_len+1.1*m;// for background studies;
  G4double ent_rin = 31.75*mm;
  G4double ent_rou = ent_rin+0.120*mm;
  
  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );
  
  G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, Mat->GetMaterial("SSteel"), "ent_log", 0, 0, 0);
  G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, Mat->GetMaterial("Vacuum"), "entvac_log", 0, 0, 0);
  
  // EFuchey: 2017/02/14
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist), entLog, "ent_phys", worldlog, false,0);
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist), entvacLog, "entvac_phys", worldlog,false,0);
     
  /*
  // EFuchey: 2017/02/14: add the possibility to change the first parameters for the beam line polycone 
  // Default set of values;
  //double z0 = sc_exbeampipeflange_dist, rin_0 = 6.20*cm, rout_0 = (6.20+0.28*2.54)*cm;
  
  int nsec = 7;
  G4double exit_z[]   = { sc_exbeampipeflange_dist, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };// -- Extended beamline for background studies (2016/09/07)

  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4double exit_rin[] = { 6.20*cm, 14.8*cm, 15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };// -- Extended beamline for background studies (2016/09/07)
  G4double exit_rou[] = { (6.20+0.28*2.54)*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };// -- Extended beamline for background studies (2016/09/07)
  
  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);
  
  G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, Mat->GetMaterial("Al"), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, Mat->GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);
  
  new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
  new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);
  
  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
  extLog->SetVisAttributes(extVisAtt);
  // extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  // extLog->SetVisAttributes(pipeVisAtt);
  */
  
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  entvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  entLog->SetVisAttributes(pipeVisAtt);
  
  //entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);
  //entLog_cut->SetVisAttributes(pipeVisAtt);
}

void NpolSBSBeamline::Place(G4LogicalVolume *worldlog) {
  
  ConstructCommonExitBeamline(worldlog);
  ConstructGMnBeamline(worldlog);
  
  //****** Place Shielding
  if( ExpType == "kNeutronExp"){

	PlaceRectangular(clamp_cone1_log,worldlog,"clamp_cone1_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(clamp_cone2_log,worldlog,"clamp_cone2_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(clamp_cone3_log,worldlog,"clamp_cone3_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(leadinmag_log,worldlog,"leadinmag_phys",0.0,0.0,(leadstart + magleadlen/2), 0.0*deg, 0.0*deg, 0.0*deg);

	G4RotationMatrix *windowshieldrm = new G4RotationMatrix();
	windowshieldrm->rotateY(-f48D48ang*pi/180);  // the inverse function needs this in radians
	G4RotationMatrix *iwindowshieldrm = new G4RotationMatrix( windowshieldrm->inverse() );
	G4ThreeVector Tm = (*iwindowshieldrm)*G4ThreeVector(19*cm, 0.0, f48D48dist-15*cm);
	PlaceRectangular(windowshield_log,worldlog,"windowshield_phys", Tm.getX(), Tm.getY(), Tm.getZ() ,0.0*deg,0.0*deg,0.0*deg);
	Tm = (*iwindowshieldrm)*G4ThreeVector(-19*cm, 0.0, f48D48dist-15*cm);
	PlaceRectangular(windowshield_log,worldlog,"windowshield_phys2", Tm.getX(), Tm.getY(), Tm.getZ(), 0.0*deg,0.0*deg,0.0*deg);
	Tm = (*iwindowshieldrm)*G4ThreeVector(0.0, gapheight/2+shieldblock3_height/2,f48D48dist-15*cm);
	PlaceRectangular(shieldblock3_log,worldlog,"windowshield_phys3", Tm.getX(), Tm.getY(), Tm.getZ(),0.0*deg,0.0*deg,0.0*deg);
	Tm = (*iwindowshieldrm)*G4ThreeVector(0.0, -gapheight/2-shieldblock3_height/2,f48D48dist-15*cm);
	PlaceRectangular(shieldblock3_log,worldlog,"windowshield_phys4", Tm.getX(), Tm.getY(), Tm.getZ(),0.0*deg,0.0*deg,0.0*deg);
  }

  if( ExpType == "kNeutronExp"){	
	PlaceRectangular(shield_cone1_log,worldlog,"shield_cone1_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
  }

  if( ExpType == "kGEp" || ExpType == "kGMn"){

	G4double inch = 2.54*cm;
	G4double xtemp = -( 5.5*inch/2.0 + 1.5*inch + (1.25/2.0+0.15)*m*tan(1.5*deg) + 2.5*cm/cos(1.5*deg) );
	G4double TargetCenter_zoffset = 6.50*inch;
	G4double z_formed_bellows = 133.2*cm - TargetCenter_zoffset;
	G4double zstart_lead_wall2 = z_formed_bellows + 76.09*inch + 1.71*inch + 15.75*inch + 1.0*inch;
	G4double zstop_lead_wall2 = 207.144*inch - TargetCenter_zoffset + 40.0*inch;
	G4double z_outer_magnetic = 182.33*cm - TargetCenter_zoffset;
	G4double zstart_lead_wall1 = z_outer_magnetic + 15*cm;
	G4double zpos_lead_wall2 = 0.5*(zstart_lead_wall2 + zstop_lead_wall2 );
	G4double xpos_lead_wall2 = -(8.0*inch + 2.5*cm + (zpos_lead_wall2 - 201.632*inch + TargetCenter_zoffset )*tan(1.5*deg));
	
	// Removed for RP-GEN PlaceRectangular(lead_wall1_log,worldlog,"lead_wall1_phys",xtemp,0.0,(zstart_lead_wall1+0.5*1.25*m),0.0*deg,1.5*deg,0.0*deg);
	PlaceRectangular(lead_wall2_log,worldlog,"lead_wall2_phys", -25*cm + xpos_lead_wall2, 0, 4.94*m /*zpos_lead_wall2*/, 0.0*deg,-0.0*deg,0.0*deg);
	
  }
  
   // End of SBS Beamline and Shielding construction.  May it rest in place.
}

