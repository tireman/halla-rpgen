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
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Cons.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "NpolMaterials.hh"
#include "NpolSBSBeamline.hh"
#include "NpolPolarimeter.hh"
#include "NpolDipole2.hh"

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
 
  if(ExpType == "kNeutronExp" || ExpType == "kGEn"){
	ConstructGEnBeamline();
	ConstructGEnLead();
	ConstructGEnClamp();
  } else if(ExpType == "kGEp"){ 
	ConstructGEpBeamline();  
	ConstructGEpLead();
  } else if(ExpType == "kGMn") {
	ConstructGEpBeamline();  
	ConstructGEpLead();
	// It appears that GMn is using similar setup as GEp for
	// beamline and shielding. W.T. (November 2018) Educated Guess
  }
  
}

NpolSBSBeamline::~NpolSBSBeamline() {}

G4String NpolSBSBeamline::GetName() {
  return G4String("SBS Beamline");
}

void NpolSBSBeamline::ConstructGEpBeamline(){

  ConstructGEnBeamline(); // all beamlines the same in G4SBS
}

void NpolSBSBeamline::ConstructGEnBeamline() {
  
  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );
  
  //We want to subtract this cylinder from the entry tube/pipe:
  G4Tubs *cut_cylinder = new G4Tubs("cut_cylinder", 0.0, swallrad, 1.0*m, 0.0*deg, 360.0*deg );
  
  G4RotationMatrix *cut_cylinder_rot = new G4RotationMatrix;
  cut_cylinder_rot->rotateX( -90.0*deg );
  
  G4SubtractionSolid *ent_tube_cut =
	new G4SubtractionSolid( "ent_tube_cut", ent_tube, cut_cylinder, cut_cylinder_rot, 
							G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );
  G4SubtractionSolid *ent_vac_cut =
	new G4SubtractionSolid( "ent_vac_cut", ent_vac, cut_cylinder, cut_cylinder_rot, 
							G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );
  
  entLog = new G4LogicalVolume(ent_tube, Mat->GetMaterial("SSteel"), "ent_log", 0, 0, 0);
  entvacLog =  new G4LogicalVolume(ent_vac, Mat->GetMaterial("HardVacuum"), "entvac_log", 0, 0, 0);
  
  entLog_cut =
	new G4LogicalVolume(ent_tube_cut, Mat->GetMaterial("SSteel"), "ent_log_cut", 0, 0, 0);
  entvacLog_cut =
	new G4LogicalVolume(ent_vac_cut, Mat->GetMaterial("HardVacuum"), "entvac_log_cut", 0, 0, 0);
  
  
  int nsec = 7;
  //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
  G4double exit_z[]   = { 162.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
  //G4double exit_z_vac[] = { 162.2*cm, 592.2*cm, 610.24*cm,610.35*cm, 1161.52*cm, 1161.53*cm,2726.46*cm };
  
  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  //G4double exit_rin[] = { 4.8*cm, 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  //G4double exit_rou[] = { 5.0*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };
  
  G4double exit_rin[] = { 6.065*2.54*cm/2., 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  G4double exit_rou[] = { (6.065/2.0+0.28)*2.54*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };
  
  
  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  //G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z_vac, exit_zero, exit_rin);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);
  
  extLog = new G4LogicalVolume(ext_cone, Mat->GetMaterial("Al"), "ext_log", 0, 0, 0);
  extvacLog = new G4LogicalVolume(ext_vac, Mat->GetMaterial("HardVacuum"), "extvac_log",0,0,0);
  
  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.1,0.9));
  extLog->SetVisAttributes(extVisAtt);
  
  extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  entvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);
  
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  
  extLog->SetVisAttributes(pipeVisAtt);
  entLog->SetVisAttributes(pipeVisAtt);
  entLog_cut->SetVisAttributes(pipeVisAtt);
  
  return;
  
}

void NpolSBSBeamline::ConstructGEnLead(){

  int nsec = 2;
  G4double clamp1_z[]   = { 162.2*cm, 228.0*cm};
  G4double clamp1_rin[] = { 5.0*cm, 10.5*cm};
  G4double clamp1_rou[] = { 25.0*cm, 25.0*cm};

  G4double clamp2_z[]   = { 2.45*m, 2.85*m,  };
  G4double clamp2_rin[] = { 11.00*cm, 12.0*cm };
  G4double clamp2_rou[] = { 25.0*cm, 25.0*cm};

  G4double clamp3_z[]   = { 4.4*m, 5.90*m,  };
  G4double clamp3_rin[] = { 16.0*cm, 17.00*cm };
  G4double clamp3_rou[] = { 25.0*cm, 25.0*cm};

  G4Polycone *clamp_cone1 = new G4Polycone("clamp_cone1", 0.0*deg, 360.0*deg, nsec, clamp1_z, clamp1_rin, clamp1_rou);
  clamp_cone1_log =
	new G4LogicalVolume( clamp_cone1, Mat->GetMaterial("Pb"), "clamp_cone1_log", 0, 0, 0 );

  G4Polycone *clamp_cone2 = new G4Polycone("clamp_cone2", 0.0*deg, 360.0*deg, nsec, clamp2_z, clamp2_rin, clamp2_rou);
  clamp_cone2_log =
	new G4LogicalVolume( clamp_cone2, Mat->GetMaterial("Pb"), "clamp_cone2_log", 0, 0, 0 );
  
  G4Polycone *clamp_cone3 = new G4Polycone("clamp_cone3", 0.0*deg, 360.0*deg, nsec, clamp3_z, clamp3_rin, clamp3_rou);
  clamp_cone3_log =
	new G4LogicalVolume( clamp_cone3, Mat->GetMaterial("Pb"), "clamp_cone3_log", 0, 0, 0 );

  // 290 -> 435 cm  box in the magnet
  G4double shield_z[]   = { 2.855*m, 4.395*m };
  G4double shield_rin[] = { 0.0, 0.0 };
  G4double shield_rou[] = { 10.50*cm, 17.*cm };

  G4Box  *leadbox = new G4Box( "leadbox",  25*cm, 15.0*cm, magleadlen/2 );
  G4Polycone *ext_cone = new G4Polycone("hollowing_tube", 0.0*deg, 360.0*deg, nsec, shield_z, shield_rin, shield_rou);

  G4SubtractionSolid *leadinmag =
	new G4SubtractionSolid("lead_w_hole", leadbox, ext_cone, 0, G4ThreeVector(0.0, 0.0, -magleadlen/2 - leadstart ) );

  leadinmag_log = new G4LogicalVolume( leadinmag, Mat->GetMaterial("Pb"), "leadinmag", 0, 0, 0 );

  ///////////  around opening of 48D48 ///////////////////////////////////////////////

  G4Box *extblocklead2 = new G4Box("extblocklead2", 10*cm, 65*cm, 10*cm );
  G4Box *shieldblock3 = new G4Box("shieldblock3", 17*cm/2, shieldblock3_height/2, 10*cm  );
  windowshield_log = new G4LogicalVolume(extblocklead2, Mat->GetMaterial("Pb"),"windowshield_log");
  shieldblock3_log = new G4LogicalVolume(shieldblock3, Mat->GetMaterial("Pb"),"shieldblock3_log");

  G4VisAttributes *leadVisAtt= new G4VisAttributes(G4Colour(0.15,0.15,0.15));
  clamp_cone1_log->SetVisAttributes(leadVisAtt);
  clamp_cone2_log->SetVisAttributes(leadVisAtt);
  clamp_cone3_log->SetVisAttributes(leadVisAtt);
  leadinmag_log->SetVisAttributes(leadVisAtt);
  windowshield_log->SetVisAttributes(leadVisAtt);
  shieldblock3_log->SetVisAttributes(leadVisAtt);

}

void NpolSBSBeamline::ConstructGEnClamp(){

  int nsec = 2;
  //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
  G4double shield_z[]   = { 2.5*m, 5.35*m };
  G4double shield_rin[] = { 8.12*cm, 14.32*cm};
  G4double shield_rou[] = { 10.11*cm, 16.33*cm };

  G4Polycone *shield_cone1 = new G4Polycone("shield_cone1", 0.0*deg, 360.0*deg, nsec, shield_z, shield_rin, shield_rou);
  shield_cone1_log =
	new G4LogicalVolume( shield_cone1, Mat->GetMaterial("Pb"), "shield_cone1_log", 0, 0, 0 );

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

  G4Box *lead_wall2 = new G4Box("lead_wall2", 5.0*cm/2.0, 1.25*m/2/*24.0*inch/2.0*/, 4.0*m/2 /*0.5*(zstop_lead_wall2 - zstart_lead_wall2)*/ );

  lead_wall2_log = new G4LogicalVolume( lead_wall2, Mat->GetMaterial("Pb"), "lead_wall2_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( 1.5*deg );

  G4cout << "Lead wall B (x,y,z) = (" << xpos_lead_wall2/cm << ", " << 0.0 << ", " << zpos_lead_wall2/cm << ")" << G4endl;
  
  lead_wall2_log->SetVisAttributes( lead_visatt );

}


void NpolSBSBeamline::Place(G4LogicalVolume *motherLV) {
  
  //****** Place Beamline *****//
  if( TargType == "kH2" || TargType == "k3He" || TargType == "kNeutTarg" ){ // modified --- W.T.
	// gas target -  1.5m in air
	PlaceRectangular(entLog,motherLV,"ent_phys",0.0, 0.0, -ent_len/2-1.5*m, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(entvacLog,motherLV,"entvac_phys",0.0, 0.0, -ent_len/2-1.5*m, 0.0*deg, 0.0*deg, 0.0*deg);
  } else {
	// Cryotarget - up against the chamber wall
	PlaceRectangular(entLog_cut,motherLV,"ent_phys",0.0, 0.0, -ent_len/2-swallrad_inner, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(entvacLog_cut,motherLV,"entvac_phys",0.0, 0.0, -ent_len/2-swallrad_inner, 0.0*deg, 0.0*deg, 0.0*deg);
  }
  PlaceRectangular(extLog,motherLV,"ext_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
  PlaceRectangular(extvacLog,motherLV,"extvac_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);

  //****** Place Shielding
  if( ExpType == "kNeutronExp"){

	PlaceRectangular(clamp_cone1_log,motherLV,"clamp_cone1_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(clamp_cone2_log,motherLV,"clamp_cone2_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(clamp_cone3_log,motherLV,"clamp_cone3_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
	PlaceRectangular(leadinmag_log,motherLV,"leadinmag_phys",0.0,0.0,(leadstart + magleadlen/2), 0.0*deg, 0.0*deg, 0.0*deg);

	G4RotationMatrix *windowshieldrm = new G4RotationMatrix();
	windowshieldrm->rotateY(-f48D48ang*pi/180);  // the inverse function needs this in radians
	G4RotationMatrix *iwindowshieldrm = new G4RotationMatrix( windowshieldrm->inverse() );
	G4ThreeVector Tm = (*iwindowshieldrm)*G4ThreeVector(19*cm, 0.0, f48D48dist-15*cm);
	PlaceRectangular(windowshield_log,motherLV,"windowshield_phys", Tm.getX(), Tm.getY(), Tm.getZ() ,0.0*deg,0.0*deg,0.0*deg);
	Tm = (*iwindowshieldrm)*G4ThreeVector(-19*cm, 0.0, f48D48dist-15*cm);
	PlaceRectangular(windowshield_log,motherLV,"windowshield_phys2", Tm.getX(), Tm.getY(), Tm.getZ(), 0.0*deg,0.0*deg,0.0*deg);
	Tm = (*iwindowshieldrm)*G4ThreeVector(0.0, gapheight/2+shieldblock3_height/2,f48D48dist-15*cm);
	PlaceRectangular(shieldblock3_log,motherLV,"windowshield_phys3", Tm.getX(), Tm.getY(), Tm.getZ(),0.0*deg,0.0*deg,0.0*deg);
	Tm = (*iwindowshieldrm)*G4ThreeVector(0.0, -gapheight/2-shieldblock3_height/2,f48D48dist-15*cm);
	PlaceRectangular(shieldblock3_log,motherLV,"windowshield_phys4", Tm.getX(), Tm.getY(), Tm.getZ(),0.0*deg,0.0*deg,0.0*deg);
  }

  if( ExpType == "kNeutronExp"){	
	PlaceRectangular(shield_cone1_log,motherLV,"shield_cone1_phys",0.0, 0.0, 0.0, 0.0*deg, 0.0*deg, 0.0*deg);
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
	
	// Removed for RP-GEN PlaceRectangular(lead_wall1_log,motherLV,"lead_wall1_phys",xtemp,0.0,(zstart_lead_wall1+0.5*1.25*m),0.0*deg,1.5*deg,0.0*deg);
	PlaceRectangular(lead_wall2_log,motherLV,"lead_wall2_phys",-0.75*m /*xpos_lead_wall2*/, 0, 5.7*m /*zpos_lead_wall2*/, 0.0*deg,1.5*deg,0.0*deg);
	
  }
  
   // End of SBS Beamline and Shielding construction.  May it rest in place.
}

