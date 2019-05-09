
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"
#include "G4ExtrudedSolid.hh"

#include "NpolSBSBeamline.hh"
#include "NpolMaterials.hh"
#include "NpolSBSTargetBuilder.hh" 

NpolSBSTargetBuilder::NpolSBSTargetBuilder(){ 

  fTargLen = 10.0*cm;
  fTargType = "kLD2";
  fTargDen = 10.5*atmosphere/(300*kelvin*k_Boltzmann);

  fTargPos = G4ThreeVector( 0, 0, 0 );
  fTargDir = G4ThreeVector( 0, 0, 1 );

  fTargDiameter = 8.0*cm; //default of 8 cm.
  
  fFlux = false;
  
  fSchamFlag = 0;
 
}

NpolSBSTargetBuilder::~NpolSBSTargetBuilder(){}

G4String NpolSBSTargetBuilder::GetName() {
  return G4String("SBS Target Builder");
}



void NpolSBSTargetBuilder::Place(G4LogicalVolume *worldlog) {

  fTargType = "kLD2"; // fDetCon->fTargType; W.T.
  fExpType = "kGMn";  // fDetCon->fExpType; W.T.

  BuildStandardScatCham( worldlog );
  
}

// EFuchey: 2017/02/10: Making a standard function to build the cryotarget itself.
// The code for building C16 and GEp are indeed almost identical.
void NpolSBSTargetBuilder::BuildStandardCryoTarget(G4LogicalVolume *motherlog, 
						 G4RotationMatrix *rot_targ, G4ThreeVector targ_offset){
  // Now let's make a cryotarget:
 
  G4double Rcell  = fTargDiameter/2.0;
  //These are assumptions. Probably should be made user-adjustable as well.
  G4double uthick = 0.1*mm;
  G4double dthick = 0.15*mm;
  G4double sthick = 0.2*mm;
  
  G4Tubs *TargetMother_solid = 
	new G4Tubs( "TargetMother_solid", 0, Rcell + sthick, (fTargLen+uthick+dthick)/2.0, 0.0, twopi );
  G4LogicalVolume *TargetMother_log = 
	new G4LogicalVolume( TargetMother_solid, NpolMaterials::GetInstance()->GetMaterial("Vacuum"), "TargetMother_log" );
  
  G4Tubs *TargetCell = new G4Tubs( "TargetCell", 0, Rcell, fTargLen/2.0, 0, twopi );
  
  G4LogicalVolume *TargetCell_log;
  
  if( fTargType == "kLH2" ){
    TargetCell_log = 
	  new G4LogicalVolume( TargetCell, NpolMaterials::GetInstance()->GetMaterial("LH2"), "TargetCell_log" );
  } else {
    TargetCell_log = 
	  new G4LogicalVolume( TargetCell, NpolMaterials::GetInstance()->GetMaterial("LD2"), "TargetCell_log" );
  }
  
  G4Tubs *TargetWall = new G4Tubs("TargetWall", Rcell, Rcell + sthick, fTargLen/2.0, 0, twopi );
  
  G4LogicalVolume *TargetWall_log = 
	new G4LogicalVolume( TargetWall, NpolMaterials::GetInstance()->GetMaterial("Al"), "TargetWall_log" );
  
  G4Tubs *UpstreamWindow = new G4Tubs("UpstreamWindow", 0, Rcell + sthick, uthick/2.0, 0, twopi );
  G4Tubs *DownstreamWindow = new G4Tubs("DownstreamWindow", 0, Rcell + sthick, dthick/2.0, 0, twopi );
  
  G4LogicalVolume *uwindow_log = 
	new G4LogicalVolume( UpstreamWindow, NpolMaterials::GetInstance()->GetMaterial("Al"), "uwindow_log" );
  G4LogicalVolume *dwindow_log = 
	new G4LogicalVolume( DownstreamWindow, NpolMaterials::GetInstance()->GetMaterial("Al"), "dwindow_log" );
  
  // Now place everything:
  // Need to fix this later: Union solid defining vacuum chamber 
  // needs to be defined with the cylinder as the first solid 
  // so that we can place the target as a daughter volume at the origin!
  
  G4double ztemp = -(fTargLen+uthick+dthick)/2.0;
  // Place upstream window:
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+uthick/2.0), uwindow_log, "uwindow_phys", TargetMother_log, false, 0 );
  // Place target and side walls:
  ztemp += uthick;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetCell_log, "TargetCell_phys", TargetMother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetWall_log, "TargetWall_phys", TargetMother_log, false, 0 );
  ztemp += fTargLen;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+dthick/2.0), dwindow_log, "dwindow_phys", TargetMother_log, false, 0 );
  
  G4double targ_zcenter = (uthick-dthick)/2.0; //position of target center relative to target mother volume
   
  //Compute position of target relative to scattering chamber:
  //The target center should be at
  //G4RotationMatrix *rot_temp = new G4RotationMatrix;
  //rot_temp = new G4RotationMatrix();
  //rot_temp->rotateX(90.0*deg);

  //for z of target center to be at zero, 
  G4double temp = targ_offset.y();
  targ_offset.setY(temp+targ_zcenter);
  
  new G4PVPlacement( rot_targ, targ_offset, TargetMother_log, "TargetMother_phys", motherlog, false, 0 );
  
  if( fFlux ){ //Make a sphere to compute particle flux:
    G4Sphere *fsph = 
	  new G4Sphere( "fsph", 1.5*fTargLen/2.0, 1.5*fTargLen/2.0+cm, 0.0*deg, 360.*deg,0.*deg, 150.*deg );
    G4LogicalVolume *fsph_log = 
	  new G4LogicalVolume( fsph, NpolMaterials::GetInstance()->GetMaterial("Air"), "fsph_log" );
    new G4PVPlacement( rot_targ, targ_offset, fsph_log, "fsph_phys", motherlog, false, 0 );
   
  }
  
  G4VisAttributes *Targ_visatt = new G4VisAttributes( G4Colour( 0.1, 0.05, 0.9 ) );
  TargetCell_log->SetVisAttributes( Targ_visatt );

  G4VisAttributes *TargWall_visatt = new G4VisAttributes( G4Colour( 0.9, .05, 0.1 ) );
  TargWall_visatt->SetForceWireframe( true );
  TargetWall_log->SetVisAttributes( TargWall_visatt );
  uwindow_log->SetVisAttributes( TargWall_visatt );
  dwindow_log->SetVisAttributes( TargWall_visatt );
  TargetMother_log->SetVisAttributes( G4VisAttributes::Invisible );
}

// ------ End of Cyro-target builder ---------- //



// ------ This function is meant to build the "Standard" scattering chamber for GMn ------- //
void NpolSBSTargetBuilder::BuildStandardScatCham(G4LogicalVolume *worldlog ){
  G4double inch = 2.54*cm;
  
  G4LogicalVolume *logicScatChamberTank =0;
  G4LogicalVolume *logicScatChamberExitFlangePlate =0;
  G4LogicalVolume *logicScatChamberFrontClamshell =0;
  G4LogicalVolume *logicScatChamberBackClamshell =0;
  G4LogicalVolume *logicScatChamberLeftSnoutWindow =0;
  G4LogicalVolume *logicScatChamberLeftSnoutWindowFrame =0;
  G4LogicalVolume *logicScatChamberRightSnoutWindow =0;
  G4LogicalVolume *logicScatChamberRightSnoutWindowFrame =0;
  G4LogicalVolume *logicScatChamber =0;

  // Scattering chamber tank:
  // basic volume:
  G4double SCHeight = 44.75*inch;
  G4double SCRadius = 20.0*inch;
  G4double SCTankThickness = 2.5*inch;
  G4double SCTankRadius = SCRadius+SCTankThickness;
  G4double SCTankHeight = SCHeight;
  G4double SCOffset = 3.75*inch;
  
  G4Tubs* solidSCTank_0 = 
    new G4Tubs("SCTank_0", SCRadius, SCTankRadius, 0.5*SCTankHeight, 0.0*deg, 360.0*deg);
  
  // exit flange:
  G4double SCExitFlangePlateHLength = 22.5*sin(25.5*atan(1)/45.0)*inch;
  G4double SCExitFlangePlateHeight = 11.0*inch;
  G4double SCExitFlangePlateThick = 1.25*inch;
  G4double SCExitFlangeHAngleApert = atan(SCExitFlangePlateHLength/(SCTankRadius+SCExitFlangePlateThick));
  G4double SCExitFlangeMaxRad = SCExitFlangePlateHLength/sin(SCExitFlangeHAngleApert);
  
  G4Tubs* solidSCExitFlangetubs = 
    new G4Tubs("SCExFlange_tubs", SCRadius, SCExitFlangeMaxRad,
			   0.5*SCExitFlangePlateHeight, 0.0, 2.0*SCExitFlangeHAngleApert);
  
  G4RotationMatrix* rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHAngleApert);
  
  G4UnionSolid* solidSCTank_0_exft = 
    new G4UnionSolid("solidSCTank_0_exft", solidSCTank_0, solidSCExitFlangetubs,
		     rot_temp, G4ThreeVector(0,0,SCOffset));
  
  G4Box* ExitFlangeHeadCut = new G4Box("ExitFlangeHeadCut", 0.5*m, 0.5*m, 0.5*m); 
    
  G4SubtractionSolid* solidSCTank_0_exf = 
    new G4SubtractionSolid("solidSCTank_0_exf", solidSCTank_0_exft, ExitFlangeHeadCut,
			   0, G4ThreeVector(-SCTankRadius-0.5*m,0,0));
  
  // exit flange hole:
  G4double SCExitFlangeHoleHeight = 7.85*inch;
  G4double SCExitFlangeHoleAngleApert = 38.25*deg;
  
  G4Tubs* solidSCExFH = 
    new G4Tubs("SCExFH", SCRadius-1.0*cm, SCTankRadius+1.5*inch,
	       0.5*SCExitFlangeHoleHeight, 0.0, SCExitFlangeHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHoleAngleApert*0.5);
  
  G4SubtractionSolid* solidSCTank_0_exfh = 
    new G4SubtractionSolid("solidSCTank_0_exfh", solidSCTank_0_exf, solidSCExFH,
			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  // windows holes: 
  G4double SCWindowHeight = 18.0*inch;
  G4double SCWindowAngleApert = 149.0*deg;
  G4double SCWindowAngleOffset = 11.0*deg;
  
  G4Tubs* solidSCWindow = 
    new G4Tubs("SCWindow", SCRadius-1.0*cm, SCTankRadius+1.0*cm, 0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wf = 
    new G4SubtractionSolid("solidSCTank_0_wf", solidSCTank_0_exfh, solidSCWindow,
			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wb = 
    new G4SubtractionSolid("solidSCTank_0_wb", solidSCTank_0_wf, solidSCWindow,
  			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  // Solid Scattering chamber tank
  G4VSolid* solidScatChamberTank = solidSCTank_0_wb;
  
  // Logic scat chamber tank
  logicScatChamberTank = 
    new G4LogicalVolume(solidScatChamberTank, NpolMaterials::GetInstance()->GetMaterial("Al"), "ScatChamberTank_log");
  
  // Scattering chamber tank placement:
  G4RotationMatrix* rotSC = new G4RotationMatrix();
  rotSC->rotateX(-90.0*deg);
  
  G4ThreeVector* SCPlacement = new G4ThreeVector(0,0,-SCOffset);
  SCPlacement->rotateX(90*deg);
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamberTank, "ScatChamberTankPhys", worldlog, false, 0);

  //Exit Flange Plate
  G4Box* solidSCExitFlangePlate = 
    new G4Box("SCExitFlangePlate_sol", SCExitFlangePlateThick/2.0,SCExitFlangePlateHeight/2.0, SCExitFlangePlateHLength); 
  
  logicScatChamberExitFlangePlate = 
    new G4LogicalVolume(solidSCExitFlangePlate, NpolMaterials::GetInstance()->GetMaterial("Al"), "SCExitFlangePlate_log");

  new G4PVPlacement(0, G4ThreeVector(-SCTankRadius-SCExitFlangePlateThick/2.0,0,0), 
  	    logicScatChamberExitFlangePlate, "SCExitFlangePlate", worldlog, false, 0); 
 
  
  // Front and back Clamshells...
  // Basic solid: 
  G4double SCClamHeight = 20.0*inch;
  G4double SCClamThick = 1.25*inch;
  G4double SCClamAngleApert = 151.0*deg;
  
  G4Tubs* solidSCClamshell_0 = 
    new G4Tubs("solidSCClamshell_0", SCTankRadius, SCTankRadius+SCClamThick, 
	       0.5*SCClamHeight, 0.0, SCClamAngleApert);
  
  // Front Clamshell:
  G4double SCFrontClamOuterRadius = SCTankRadius+SCClamThick;
  G4double SCBeamExitAngleOffset = 64.5*deg;
  G4double SCLeftSnoutAngle = -24.2*deg;
  G4double SCLeftSnoutAngleOffset = SCBeamExitAngleOffset+SCLeftSnoutAngle;
  G4double SCRightSnoutAngle = 50.1*deg;
  G4double SCRightSnoutAngleOffset = SCBeamExitAngleOffset+SCRightSnoutAngle;
  
  // Snouts: NB: similarly to the GEp scattering chamber, 
  // the "left" and "right" is defined as viewed from downstream.
  // In other words, the "left" snout is actually on the right side of the beam 
  // (looking on the right direction) and vice-versa.
  
  // Right snout opening:
  G4double SCRightSnoutDepth = 15.0*inch;// x
  G4double SCRightSnoutWidth = 26.0*inch;// y
  G4double SCRightSnoutHeight = 18.0*inch; //z
  G4double SCRightSnoutHoleHeight = 14.0*inch;
  G4double SCRightSnoutHoleAngleApert = 50.0*deg;
  G4double SCRightSnoutBoxAngleOffset = 9.40*deg;
  
  // Basic "box"
  G4Box* RightSnoutBox = 
    new G4Box("RightSnoutBox", SCRightSnoutDepth*0.5, SCRightSnoutWidth*0.5, SCRightSnoutHeight*0.5);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg-SCRightSnoutAngleOffset-SCRightSnoutBoxAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_rsb = 
    new G4UnionSolid("solidSCClamshell_0_rsb", solidSCClamshell_0, 
		     RightSnoutBox, rot_temp, 
		     G4ThreeVector(SCFrontClamOuterRadius*sin(90.0*deg-SCRightSnoutAngleOffset),
		      		   SCFrontClamOuterRadius*cos(90.0*deg-SCRightSnoutAngleOffset), 
		     		   0) );
  
  // Basic box cut: remove the surplus outside of the clamshell
  // NB: the surplus inside is removed along with the one of the left snout
  G4Box* RightSnoutBox_cut = 
    new G4Box("RightSnoutBox_cut", 
	      SCRightSnoutDepth+0.01*inch, SCRightSnoutWidth, SCRightSnoutHeight*0.5+0.01*inch);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg-SCRightSnoutAngleOffset);
  
  G4SubtractionSolid* solidSCFrontClam_0_rsbc = 
    new G4SubtractionSolid("solidSCFrontClam_0_rsbc", solidSCFrontClam_0_rsb,
			   RightSnoutBox_cut, rot_temp, 
			   G4ThreeVector((SCFrontClamOuterRadius+SCRightSnoutDepth)
					 *sin(90.0*deg-SCRightSnoutAngleOffset), 
					 (SCFrontClamOuterRadius+SCRightSnoutDepth)
					 *cos(90.0*deg-SCRightSnoutAngleOffset),
					 0) );
  
  // Cutting the hole...
  G4Tubs* RightSnoutApertCut = 
    new G4Tubs("RightSnoutApertCut", 0.0, 30.0*inch, SCRightSnoutHoleHeight/2.0, 
	       0.0, SCRightSnoutHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-SCRightSnoutAngleOffset+SCRightSnoutHoleAngleApert*0.5);
  
  G4SubtractionSolid* solidSCFrontClam_0_rs =
    new G4SubtractionSolid("solidSCClamshell_0_rs", solidSCFrontClam_0_rsbc, 
			   RightSnoutApertCut, rot_temp, G4ThreeVector(0, 0, 0));
  
  
  // Right snout window+frame:
  G4double SCRightSnoutWindowWidth = 26.351*inch;
  G4double SCSnoutWindowThick = 0.02*inch;
  G4double SCSnoutWindowFrameThick = 0.75*inch;
  
  G4double SCRightSnoutWindowDist = 23.74*inch+SCSnoutWindowThick*0.5;
  G4double SCRightSnoutHoleWidth = 21.855*inch;
  G4double SCRightSnoutHoleCurvRad = 2.1*inch;
  G4double SCRightSnoutWindowFrameDist = SCRightSnoutWindowDist+SCSnoutWindowThick*0.5+SCSnoutWindowFrameThick*0.5;
  
  // window
  G4Box* solidRightSnoutWindow = 
    new G4Box("RightSnoutWindow_sol", SCRightSnoutWindowWidth*0.5, 
	      SCRightSnoutHeight*0.5, SCSnoutWindowThick*0.5);

  logicScatChamberRightSnoutWindow = 
    new G4LogicalVolume(solidRightSnoutWindow, NpolMaterials::GetInstance()->GetMaterial("Al"), "SCRightSnoutWindow_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCRightSnoutAngle);
  
  new G4PVPlacement(rot_temp, 
  	    G4ThreeVector(SCRightSnoutWindowDist*sin(SCRightSnoutAngle),0,
					  SCRightSnoutWindowDist*cos(SCRightSnoutAngle)),
					logicScatChamberRightSnoutWindow, "SCRightSnoutWindow", worldlog, false, 0);
    
  // basic window frame
  G4Box* solidRightSnoutWindowFrame_0 = 
    new G4Box("RightSnoutWindowFrame_0", SCRightSnoutWidth*0.5, 
	      SCRightSnoutHeight*0.5, SCSnoutWindowFrameThick*0.5);
  
  // + lots of cut out solids...
  G4Tubs* solidRightSnoutWindowFrame_cut_0 = 
    new G4Tubs("RightSnoutWindowFrame_cut_0", 0.0, SCRightSnoutHoleCurvRad, SCSnoutWindowFrameThick, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_0_htr = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_htr", solidRightSnoutWindowFrame_0,
			   solidRightSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(SCRightSnoutHoleWidth*0.5-SCRightSnoutHoleCurvRad, 
					 SCRightSnoutHoleHeight*0.5-SCRightSnoutHoleCurvRad,
					 0) );

  G4SubtractionSolid* solidRightSnoutWindowFrame_0_hbr = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_hbr", solidRightSnoutWindowFrame_0_htr,
			   solidRightSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(SCRightSnoutHoleWidth*0.5-SCRightSnoutHoleCurvRad, 
					 -SCRightSnoutHoleHeight*0.5+SCRightSnoutHoleCurvRad,
					 0) );
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_0_hbl = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_hbl", solidRightSnoutWindowFrame_0_hbr,
			   solidRightSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(-SCRightSnoutHoleWidth*0.5+SCRightSnoutHoleCurvRad, 
					 -SCRightSnoutHoleHeight*0.5+SCRightSnoutHoleCurvRad,
					 0) );
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_0_htl = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame_0_htl", solidRightSnoutWindowFrame_0_hbl,
			   solidRightSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(-SCRightSnoutHoleWidth*0.5+SCRightSnoutHoleCurvRad, 
					 SCRightSnoutHoleHeight*0.5-SCRightSnoutHoleCurvRad,
					 0));
  
  G4Box* solidRightSnoutWindowFrame_cut_1 = 
    new G4Box("RightSnoutWindowFrame_cut_1", SCRightSnoutHoleWidth*0.5-SCRightSnoutHoleCurvRad, 
	      SCRightSnoutHoleHeight*0.5, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidRightSnoutWindowFrame_1 = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame", solidRightSnoutWindowFrame_0_htl,
						   solidRightSnoutWindowFrame_cut_1, 0, G4ThreeVector(0,0,0));
  
  G4Box* solidRightSnoutWindowFrame_cut_2 = 
    new G4Box("RightSnoutWindowFrame_cut_2", SCRightSnoutHoleWidth*0.5,
	      SCRightSnoutHoleHeight*0.5-SCRightSnoutHoleCurvRad, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidRightSnoutWindowFrame = 
    new G4SubtractionSolid("solidRightSnoutWindowFrame", solidRightSnoutWindowFrame_1,
  			   solidRightSnoutWindowFrame_cut_2, 0, G4ThreeVector(0,0,0));

  logicScatChamberRightSnoutWindowFrame = 
    new G4LogicalVolume(solidRightSnoutWindowFrame, 
						NpolMaterials::GetInstance()->GetMaterial("Al"), "SCRightSnoutWindowFrame_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCRightSnoutAngle);

  /*new G4PVPlacement(rot_temp,G4ThreeVector(SCRightSnoutWindowFrameDist*sin(SCRightSnoutAngle),0,
  	  SCRightSnoutWindowFrameDist*cos(SCRightSnoutAngle)), 
	  logicScatChamberRightSnoutWindowFrame, "SCRightSnoutWindowFrame", worldlog, false, 0);*/
  // W.T.  Dying here!!! 
		    
  // Left snout opening:
  G4double SCLeftSnoutDepth = 4.0*inch;// x
  G4double SCLeftSnoutWidth = 16.338*inch;// y
  G4double SCLeftSnoutHeight = 11.0*inch; //z
  G4double SCLeftSnoutHoleHeight = 7.0*inch;
  G4double SCLeftSnoutHoleAngleApert = 30.0*deg;
  G4double SCLeftSnoutYOffset = 0.50*inch;
  
  // Basic box
  G4Box* LeftSnoutBox = 
    new G4Box("LeftSnoutBox", SCLeftSnoutDepth*0.5, SCLeftSnoutWidth*0.5, SCLeftSnoutHeight*0.5);

  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg-SCLeftSnoutAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_lsb = 
    new G4UnionSolid("solidSCClamshell_0_lsb", solidSCFrontClam_0_rs, 
		     LeftSnoutBox, rot_temp, 
		     G4ThreeVector((SCFrontClamOuterRadius-1.85*inch)*sin(90.0*deg-SCLeftSnoutAngleOffset),
		     		   (SCFrontClamOuterRadius-1.85*inch)*cos(90.0*deg-SCLeftSnoutAngleOffset), 
		     		   SCLeftSnoutYOffset) );
  
  // remove all surplus material inside the scat chamber
  G4Tubs* SnoutsInnerCut = 
    new G4Tubs("SnoutsInnerCut", 0.0, SCTankRadius, SCClamHeight/2.0, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCFrontClam_0_lsbc =
    new G4SubtractionSolid("solidSCClamshell_0_lsbc", solidSCFrontClam_0_lsb, 
			   SnoutsInnerCut, 0, G4ThreeVector());
  
  // Cut the hole
  G4Tubs* LeftSnoutApertCut = 
    new G4Tubs("LeftSnoutApertCut", 0.0, 30.0*inch, SCLeftSnoutHoleHeight/2.0, 
	       0.0, SCLeftSnoutHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-SCLeftSnoutAngleOffset+SCLeftSnoutHoleAngleApert*0.5);
 
  G4SubtractionSolid* solidSCFrontClam_0_ls =
    new G4SubtractionSolid("solidSCClamshell_0_ls", solidSCFrontClam_0_lsbc, 
			   LeftSnoutApertCut, rot_temp, G4ThreeVector(0, 0, SCLeftSnoutYOffset));
  
  // Left snout window+frame:
  G4double SCLeftSnoutWindowDist = 23.9*inch+SCSnoutWindowThick*0.5;
  G4double SCLeftSnoutWindowFrameDist = SCLeftSnoutWindowDist+SCSnoutWindowThick*0.5+SCSnoutWindowFrameThick*0.5;
  G4double SCLeftSnoutHoleWidth = 12.673*inch;
  G4double SCLeftSnoutHoleCurvRad = 1.05*inch;
  
  // window
  G4Box* solidLeftSnoutWindow_0 = 
    new G4Box("LeftSnoutWindow_sol", SCLeftSnoutWidth*0.5, SCLeftSnoutHeight*0.5, SCSnoutWindowThick*0.5);
  
  G4Tubs* solidLeftSnoutWindow_cut_0 = 
    new G4Tubs("LeftSnoutWindow_cut_0", 0.0, 2.579*inch, SCSnoutWindowFrameThick, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidLeftSnoutWindow = 
    new G4SubtractionSolid("solidLeftSnoutWindow", solidLeftSnoutWindow_0, solidLeftSnoutWindow_cut_0,
			   0, G4ThreeVector(10.287*inch, -SCLeftSnoutYOffset, 0));
  
  logicScatChamberLeftSnoutWindow = 
    new G4LogicalVolume(solidLeftSnoutWindow, NpolMaterials::GetInstance()->GetMaterial("Al"), "SCLeftSnoutWindow_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCLeftSnoutAngle);
  
  new G4PVPlacement(rot_temp,G4ThreeVector(SCLeftSnoutWindowDist*sin(SCLeftSnoutAngle),
  			  SCLeftSnoutYOffset, SCLeftSnoutWindowDist*cos(SCLeftSnoutAngle)), 
  	    logicScatChamberLeftSnoutWindow, "SCLeftSnoutWindow", worldlog, false, 0);
 
  // window frame
  G4Box* solidLeftSnoutWindowFrame_0 = 
    new G4Box("LeftSnoutWindowFrame_0", SCLeftSnoutWidth*0.5, 
	      SCLeftSnoutHeight*0.5, SCSnoutWindowFrameThick*0.5);
  
  // + lots of cut out solids...
  G4Tubs* solidLeftSnoutWindowFrame_cut_0 = 
    new G4Tubs("LeftSnoutWindowFrame_cut_0", 0.0, SCLeftSnoutHoleCurvRad, SCSnoutWindowFrameThick, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_htr = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_htr", solidLeftSnoutWindowFrame_0,
			   solidLeftSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(SCLeftSnoutHoleWidth*0.5-SCLeftSnoutHoleCurvRad, 
					 SCLeftSnoutHoleHeight*0.5-SCLeftSnoutHoleCurvRad,
					 0) );

  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_hbr = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_hbr", solidLeftSnoutWindowFrame_0_htr,
			   solidLeftSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(SCLeftSnoutHoleWidth*0.5-SCLeftSnoutHoleCurvRad, 
					 -SCLeftSnoutHoleHeight*0.5+SCLeftSnoutHoleCurvRad,
					 0) );
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_hbl = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_hbl", solidLeftSnoutWindowFrame_0_hbr,
			   solidLeftSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(-SCLeftSnoutHoleWidth*0.5+SCLeftSnoutHoleCurvRad, 
					 -SCLeftSnoutHoleHeight*0.5+SCLeftSnoutHoleCurvRad,
					 0) );
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_0_htl = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_0_htl", solidLeftSnoutWindowFrame_0_hbl,
			   solidLeftSnoutWindowFrame_cut_0, 0, 
			   G4ThreeVector(-SCLeftSnoutHoleWidth*0.5+SCLeftSnoutHoleCurvRad, 
					 SCLeftSnoutHoleHeight*0.5-SCLeftSnoutHoleCurvRad,
					 0) );
  
  G4Box* solidLeftSnoutWindowFrame_cut_1 = 
    new G4Box("LeftSnoutWindowFrame_cut_1", SCLeftSnoutHoleWidth*0.5-SCLeftSnoutHoleCurvRad, 
	      SCLeftSnoutHoleHeight*0.5, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_1 = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame", solidLeftSnoutWindowFrame_0_htl,
  			   solidLeftSnoutWindowFrame_cut_1, 0, G4ThreeVector());
  
  G4Box* solidLeftSnoutWindowFrame_cut_2 = 
    new G4Box("LeftSnoutWindowFrame_cut_2", SCLeftSnoutHoleWidth*0.5,
	      SCLeftSnoutHoleHeight*0.5-SCLeftSnoutHoleCurvRad, SCSnoutWindowFrameThick);
  
  G4SubtractionSolid* solidLeftSnoutWindowFrame_2 = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame_2", solidLeftSnoutWindowFrame_1,
  			   solidLeftSnoutWindowFrame_cut_2, 0, G4ThreeVector());

  G4SubtractionSolid* solidLeftSnoutWindowFrame = 
    new G4SubtractionSolid("solidLeftSnoutWindowFrame", solidLeftSnoutWindowFrame_2, solidLeftSnoutWindow_cut_0,
			   0, G4ThreeVector(10.287*inch, -SCLeftSnoutYOffset, 0));
  
  logicScatChamberLeftSnoutWindowFrame = 
    new G4LogicalVolume(solidLeftSnoutWindowFrame, 
						NpolMaterials::GetInstance()->GetMaterial("Al"), "SCLeftSnoutWindowFrame_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateY(-SCLeftSnoutAngle);
  
  new G4PVPlacement(rot_temp, G4ThreeVector(SCLeftSnoutWindowFrameDist*sin(SCLeftSnoutAngle),
  			  SCLeftSnoutYOffset, SCLeftSnoutWindowFrameDist*cos(SCLeftSnoutAngle)), 
  	    logicScatChamberLeftSnoutWindowFrame, "SCLeftSnoutWindow", worldlog, false, 0);

  // Exit Beam Pipe:
  // Should come after left snout
  G4Tubs* solidExitBeamPipe = new G4Tubs("solidExitBeamPipe", 
					 0.0, 55.0*mm, 2.150*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg+SCBeamExitAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_ebp =
    new G4UnionSolid("solidSCClamshell_0_ebp", solidSCFrontClam_0_ls, 
		     solidExitBeamPipe, rot_temp, 
		     G4ThreeVector((SCFrontClamOuterRadius+2.003*inch)*sin(90.0*deg-SCBeamExitAngleOffset),
				   (SCFrontClamOuterRadius+2.003*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
				   0.0) );
  
  // remove extra material from the left snout around the chamber
  G4Tubs* solidExitBeamPipeSurroundCut = new G4Tubs("solidExitBeamPipeSurroundCut", 
						    55.0*mm, 2.803*inch, 1.684*inch, 0.0, 360.0*deg);
  
   G4SubtractionSolid* solidSCFrontClam_0_ebps =
    new G4SubtractionSolid("solidSCClamshell_0_ebps", solidSCFrontClam_0_ebp, 
		     solidExitBeamPipeSurroundCut, rot_temp, 
		     G4ThreeVector((SCFrontClamOuterRadius+1.684*inch)*sin(90.0*deg-SCBeamExitAngleOffset),
				   (SCFrontClamOuterRadius+1.684*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
				   0.0) );
   
  G4Tubs* solidExitBeamPipeFlange = new G4Tubs("solidExitBeamPipeFlange", 
					       0.0, 2.985*inch, 0.3925*inch, 0.0, 360.0*deg);
  
  G4UnionSolid* solidSCFrontClam_0_ebpf =
    new G4UnionSolid("solidSCClamshell_0_ebpf", solidSCFrontClam_0_ebps, 
		     solidExitBeamPipeFlange, rot_temp, 
		     G4ThreeVector((SCFrontClamOuterRadius+3.7605*inch)*sin(90.0*deg-SCBeamExitAngleOffset), 
				   (SCFrontClamOuterRadius+3.7605*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
				   0.0) );
  
  G4Tubs* solidExitBeamPipeHole = new G4Tubs("solidBackViewPipeHole", 
					     0.0, 50.0*mm, 7.903*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCFrontClam_0_ebph =
    new G4SubtractionSolid("solidSCClamshell_0_ebph", solidSCFrontClam_0_ebpf, 
			   solidExitBeamPipeHole, rot_temp, 
			   G4ThreeVector(SCFrontClamOuterRadius*sin(90.0*deg-SCBeamExitAngleOffset),
					 SCFrontClamOuterRadius*cos(90.0*deg-SCBeamExitAngleOffset), 
					 0.0) );
  
  // Placing Front ClamShell (at last)
  G4VSolid* solidSCFrontClamShell = solidSCFrontClam_0_ebph;
  
  logicScatChamberFrontClamshell = 
    new G4LogicalVolume(solidSCFrontClamShell, NpolMaterials::GetInstance()->GetMaterial("Al"), "SCFrontClamshell_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateZ(90.0*deg+SCClamAngleApert*0.5-SCWindowAngleOffset);
  
  new G4PVPlacement(rot_temp,G4ThreeVector(0,0,0),logicScatChamberFrontClamshell,"SCFrontClamshell", worldlog, false, 0);
  
  // Back Clamshell:
  G4double SCBackClamThick = 0.80*inch;
  G4double SCBackClamHeight = 12.0*inch;
  G4double SCBackClamAngleApert = 145.7*deg;
  G4double SCBackClamAngleOffset = -2.5*deg;
  G4double SCBackClamOuterRadius = SCTankRadius+SCClamThick+SCBackClamThick;
  G4double SCBeamEntranceAngleOffset = SCBackClamAngleOffset-84*deg;
  G4double SCBackViewPipeAngleOffset = SCBackClamAngleOffset-39.9*deg;
  
  G4Tubs* solidSCAddBackClam = 
    new G4Tubs("solidSCAddBackClam", SCTankRadius+SCClamThick-0.5*cm, SCBackClamOuterRadius, 
	       0.5*SCBackClamHeight, 0.0, SCBackClamAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(SCBackClamAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_abc = 
    new G4UnionSolid("solidSCClamshell_0_abc", solidSCClamshell_0, solidSCAddBackClam,
		     rot_temp, G4ThreeVector(0,0,0));
  
  // Entrance beam pipe
  G4Tubs* solidEntranceBeamPipe = 
    new G4Tubs("solidEntranceBeamPipe", 0.0, 2.25*inch, 0.825*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg-SCBeamEntranceAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_ebp =
    new G4UnionSolid("solidSCClamshell_0_ebp", solidSCBackClam_0_abc, solidEntranceBeamPipe, rot_temp, 
		     G4ThreeVector(SCBackClamOuterRadius*sin(90.0*deg+SCBeamEntranceAngleOffset),
				   SCBackClamOuterRadius*cos(90.0*deg+SCBeamEntranceAngleOffset),
				   0.0) );
  
  G4Tubs* solidEntranceBeamPipeHole = 
    new G4Tubs("solidEntranceBeamPipeHole", 0.0, 1.0*inch, 5.375*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCBackClam_0_ebph =
    new G4SubtractionSolid("solidSCClamshell_0_ebph", solidSCBackClam_0_ebp, 
			   solidEntranceBeamPipeHole, rot_temp, 
			   G4ThreeVector(SCBackClamOuterRadius*
					 sin(90.0*deg+SCBeamEntranceAngleOffset),
					 SCBackClamOuterRadius*
					 cos(90.0*deg+SCBeamEntranceAngleOffset),
					 0.0) );
  
  // Back view pipe
  G4Tubs* solidBackViewPipe = new G4Tubs("solidBackViewPipe", 
					 0.0, 2.12*inch, 1.555*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg-SCBackViewPipeAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_bvp =
    new G4UnionSolid("solidSCClamshell_0_bvp", solidSCBackClam_0_ebph, 
		     solidBackViewPipe, rot_temp, 
		     G4ThreeVector((SCBackClamOuterRadius+1.501*inch)*sin(90.0*deg+SCBackViewPipeAngleOffset),
				   (SCBackClamOuterRadius+1.501*inch)*cos(90.0*deg+SCBackViewPipeAngleOffset), 
				   0.0) );
  
  G4Tubs* solidBackViewPipeFlange = new G4Tubs("solidBackViewPipeFlange", 
					       0.0, 3.0*inch, 0.5*inch, 0.0, 360.0*deg);
  
  G4UnionSolid* solidSCBackClam_0_bvpf =
    new G4UnionSolid("solidSCClamshell_0_bvpf", solidSCBackClam_0_bvp, 
		     solidBackViewPipeFlange, rot_temp, 
		     G4ThreeVector((SCBackClamOuterRadius+2.556*inch)*sin(90.0*deg+SCBackViewPipeAngleOffset), 
				   (SCBackClamOuterRadius+2.556*inch)*cos(90.0*deg+SCBackViewPipeAngleOffset), 
				   0.0) );
  
  G4Tubs* solidBackViewPipeHole = new G4Tubs("solidBackViewPipeHole", 
					     0.0, 2.0*inch, 4.0*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCBackClam_0_bvph =
    new G4SubtractionSolid("solidSCClamshell_0_bvph", solidSCBackClam_0_bvpf, 
			   solidBackViewPipeHole, rot_temp, 
			   G4ThreeVector(SCBackClamOuterRadius*sin(90.0*deg+SCBackViewPipeAngleOffset),
					 SCBackClamOuterRadius*cos(90.0*deg+SCBackViewPipeAngleOffset), 
					 0.0) );
  
  
  // Placing back clamshell
  G4VSolid* solidSCBackClamShell = solidSCBackClam_0_bvph;
    
  logicScatChamberBackClamshell = new G4LogicalVolume(solidSCBackClamShell, 
						      NpolMaterials::GetInstance()->GetMaterial("Al"), 
						      "SCBackClamshell_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateZ(-90.0*deg+SCClamAngleApert*0.5+SCWindowAngleOffset);
  
  new G4PVPlacement(rot_temp,G4ThreeVector(0,0,0),logicScatChamberBackClamshell,"SCBackClamshell", worldlog, false, 0);
  
  // Scattering chamber volume
  //
  G4Tubs* solidScatChamber_0 = new G4Tubs("SC", 0.0, SCRadius, 0.5* SCHeight, 0.0*deg, 360.0*deg);
  
  G4Tubs* solidSCWindowVacuum = 
    new G4Tubs("SCWindowVacuumFront", SCRadius-1.0*cm, SCTankRadius, 0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4UnionSolid* solidScatChamber_0_wbv = 
    new G4UnionSolid("solidScatChamber_0_wbv", solidScatChamber_0, solidSCWindowVacuum,
			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  //G4UnionSolid* solidScatChamber_0_wfv = 
  //new G4UnionSolid("solidScatChamber_0_wfv", solidScatChamber_0_wbv, solidSCWindowVacuum,
  //		   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  G4UnionSolid* solidScatChamber_0_entbp = 
    new G4UnionSolid("solidScatChamber_0_entbp", solidScatChamber_0_wbv, solidEntranceBeamPipeHole, 
		     rot_temp, G4ThreeVector(0, -SCRadius, SCOffset));
  //solidExitBeamPipeHole
  
  G4UnionSolid* solidScatChamber_0_exbp = 
    new G4UnionSolid("solidScatChamber_0_exbp", solidScatChamber_0_entbp, solidExitBeamPipeHole, 
		     rot_temp, G4ThreeVector(0, +SCRadius, SCOffset));
  
  G4VSolid* solidScatChamber = solidScatChamber_0_exbp;

  logicScatChamber = new G4LogicalVolume(solidScatChamber, NpolMaterials::GetInstance()->GetMaterial("Vacuum"), "ScatChamber_log");
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamber, "ScatChamberPhys", worldlog, false, 0);
 
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  
  //Call BuildStandardCryoTarget HERE !
  BuildStandardCryoTarget(logicScatChamber, rot_temp, G4ThreeVector(0, 0, SCOffset));
  
  G4VisAttributes* Invisible  = new G4VisAttributes(G4Colour(0.,0.,0.)); 
  Invisible->SetVisibility(false);
  G4VisAttributes* colourDarkGrey = new G4VisAttributes(G4Colour(0.3,0.3,0.3)); 
  colourDarkGrey->SetForceWireframe(true);
  G4VisAttributes* colourGrey = new G4VisAttributes(G4Colour(0.7,0.7,0.7)); 
  colourGrey->SetForceWireframe(true);
  G4VisAttributes* colourCyan = new G4VisAttributes(G4Colour(0.,1.,1.)); 
  
  logicScatChamberTank->SetVisAttributes(colourDarkGrey);
  logicScatChamberFrontClamshell->SetVisAttributes(colourGrey);
  logicScatChamberLeftSnoutWindow->SetVisAttributes(Invisible);
  logicScatChamberLeftSnoutWindowFrame->SetVisAttributes(colourCyan);
  logicScatChamberRightSnoutWindow->SetVisAttributes(Invisible);
  logicScatChamberRightSnoutWindowFrame->SetVisAttributes(colourCyan);
  logicScatChamberBackClamshell->SetVisAttributes(colourGrey);
  logicScatChamberExitFlangePlate->SetVisAttributes(colourGrey);
  logicScatChamber->SetVisAttributes(Invisible);
 
}


void NpolSBSTargetBuilder::BuildGasTarget(G4LogicalVolume *worldlog){

  if( fFlux ){ //Make a sphere to compute particle flux:
    G4Sphere *fsph = new G4Sphere( "fsph", 1.5*fTargLen/2.0, 1.5*fTargLen/2.0+cm, 0.0*deg, 360.*deg,
				   0.*deg, 150.*deg );
    G4LogicalVolume *fsph_log = new G4LogicalVolume( fsph, NpolMaterials::GetInstance()->GetMaterial("Air"), "fsph_log" );
    new G4PVPlacement( 0, G4ThreeVector(0,0,0), fsph_log, "fsph_phys", worldlog, false, 0 );

    fsph_log->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    
    /*G4String FluxSDname = "FLUX";
    G4String Fluxcollname = "FLUXHitsCollection";
    NpolSBSCalSD *FluxSD = NULL;
    if( !( FluxSD = (NpolSBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(FluxSDname) ) ){
      G4cout << "Adding FLUX SD to SDman..." << G4endl;
      FluxSD = new NpolSBSCalSD( FluxSDname, Fluxcollname );
      fDetCon->fSDman->AddNewDetector( FluxSD );
      (fDetCon->SDlist).insert( FluxSDname );
      fDetCon->SDtype[FluxSDname] = kCAL;

      (FluxSD->detmap).depth = 0;
    }
    fsph_log->SetSensitiveDetector( FluxSD );*/
  }
  
  //Desired minimum scattering angle for SBS = 5 deg (for SIDIS @10 deg):
  //double sbs_scattering_angle_min = 5.0*deg;
  double sc_exitpipe_radius_inner = 48.0*mm;
  double sc_exitpipe_radius_outer = 50.0*mm;
  double sc_entrypipe_radius_inner = 31.75*mm;
  //double sc_entrypipe_radius_outer = sc_entrypipe_radius_inner+0.12*mm;
  double sc_winthick = 0.38*mm;

  double zstart_sc = -1.5*m;
  double zend_sc   = 162.2*cm;

  double zpos_sc = (zstart_sc + zend_sc)/2.0;
  double dz_sc = (zend_sc - zstart_sc)/2.0 - sc_winthick;

  //Material definition was moved to ConstructMaterials();
  //--------- Glass target cell -------------------------------

  double wallthick = 1.61*mm;
  double capthick  = 0.126*mm;
  
  double sc_radius_inner = sc_exitpipe_radius_outer;
  double sc_radius_outer = sc_radius_inner + sc_winthick;

  G4Tubs *sc_wall = new G4Tubs("sc_tube", sc_radius_inner, sc_radius_outer, dz_sc, 0.*deg, 360.*deg );
  G4Tubs *sc_vacuum = new G4Tubs("sc_vac", 0.0, sc_radius_inner, dz_sc, 0.*deg, 360.*deg );

  G4Tubs *sc_cap_upstream = new G4Tubs("sc_cap_upstream", sc_entrypipe_radius_inner, sc_radius_outer, sc_winthick/2.0, 0.*deg, 360.*deg );
  G4Tubs *sc_cap_downstream = new G4Tubs("sc_cap_downstream", sc_exitpipe_radius_inner, sc_radius_outer, sc_winthick/2.0, 0.*deg, 360.*deg );
  
  G4LogicalVolume *sc_wall_log = new G4LogicalVolume( sc_wall, NpolMaterials::GetInstance()->GetMaterial("Al"), "sc_wall_log" );
  G4LogicalVolume *sc_vacuum_log = new G4LogicalVolume( sc_vacuum, NpolMaterials::GetInstance()->GetMaterial("Vacuum"), "sc_vacuum_log" );
  G4LogicalVolume *sc_cap_upstream_log = new G4LogicalVolume( sc_cap_upstream, NpolMaterials::GetInstance()->GetMaterial("Al"), "sc_cap_upstream_log" );
  G4LogicalVolume *sc_cap_downstream_log = new G4LogicalVolume( sc_cap_downstream, NpolMaterials::GetInstance()->GetMaterial("Al"), "sc_cap_downstream_log" );
  
  G4Tubs *targ_tube = new G4Tubs("targ_tube", fTargDiameter/2.0-wallthick, fTargDiameter/2.0, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, fTargDiameter/2.0, capthick/2.0, 0.*deg, 360.*deg );

  G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, NpolMaterials::GetInstance()->GetMaterial("GE180"),"targ_tube_log");
  G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, NpolMaterials::GetInstance()->GetMaterial("GE180"),"targ_cap_log");

  // gas
  G4Tubs *gas_tube = new G4Tubs("gas_tube", 0.0, fTargDiameter/2.0-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* gas_tube_log = NULL;


  if( fTargType == "kH2" || fTargType == "kNeutTarg" ){
    gas_tube_log = new G4LogicalVolume(gas_tube, NpolMaterials::GetInstance()->GetMaterial("refH2"), "gas_tube_log");
  }
  if( fTargType == "kD2" ){
    gas_tube_log = new G4LogicalVolume(gas_tube, NpolMaterials::GetInstance()->GetMaterial("refD2"), "gas_tube_log");
  }
  if( fTargType == "k3He" ){
    gas_tube_log = new G4LogicalVolume(gas_tube, NpolMaterials::GetInstance()->GetMaterial("pol3He"), "gas_tube_log");
  }

  G4LogicalVolume *motherlog = worldlog;
  double target_zpos = 0.0;
  
  if( fSchamFlag == 1 ){
    motherlog = sc_vacuum_log;
    target_zpos = -zpos_sc;
  }

  //if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos), targ_tube_log,
		    "targ_tube_phys", motherlog, false, 0);
  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos+fTargLen/2.0+capthick/2.0), targ_cap_log,
		    "targ_cap_phys1", motherlog, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos-fTargLen/2.0-capthick/2.0), targ_cap_log,
		    "targ_cap_phys2", motherlog, false, 1);
  
  assert(gas_tube_log);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos), gas_tube_log,
		    "gas_tube_phys", motherlog, false, 0);
  
  //Place scattering chamber:
  if( fSchamFlag == 1 ){
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc), sc_wall_log, "sc_wall_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc-dz_sc-sc_winthick/2.0), sc_cap_upstream_log, "sc_cap_upstream_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc+dz_sc+sc_winthick/2.0), sc_cap_downstream_log, "sc_cap_downstream_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc), sc_vacuum_log, "sc_vacuum_phys", worldlog, false, 0 );
  }
  
  sc_vacuum_log->SetVisAttributes( G4VisAttributes::Invisible );

  //Visualization attributes:
  //ScatteringChamber_log->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *sc_wall_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  sc_wall_visatt->SetForceWireframe(true);

  sc_wall_log->SetVisAttributes( sc_wall_visatt );
  sc_cap_upstream_log->SetVisAttributes( sc_wall_visatt );
  sc_cap_upstream_log->SetVisAttributes( sc_wall_visatt );

  // G4VisAttributes *sc_exit_pipe_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  // sc_exit_pipe_log->SetVisAttributes( sc_exit_pipe_visatt );

  // G4VisAttributes *sc_win_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.0 ) );
  // sc_window_sbs_log->SetVisAttributes( sc_win_visatt );
  // //sc_window_bb_log->SetVisAttributes( sc_win_visatt );

  G4VisAttributes *tgt_cell_visatt = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
  tgt_cell_visatt->SetForceWireframe(true);

  targ_cap_log->SetVisAttributes( tgt_cell_visatt );
  targ_tube_log->SetVisAttributes( tgt_cell_visatt );

  G4VisAttributes *tgt_gas_visatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 1.0 ) );
  gas_tube_log->SetVisAttributes( tgt_gas_visatt );

}
