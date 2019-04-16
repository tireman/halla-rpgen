

#include "NpolEventPreProcessing.hh"

double NpolEventPreProcessing::NpolAng = 0.431096; /*radians; angle of NPOL relative to beam axis*/

NpolEventPreProcessing *NpolEventPreProcessing::PreProcess = NULL;

NpolEventPreProcessing *NpolEventPreProcessing::GetInstance() {
	if(PreProcess == NULL)
		PreProcess = new NpolEventPreProcessing();

	return PreProcess;
}

NpolEventPreProcessing::NpolEventPreProcessing(){}

NpolEventPreProcessing::~NpolEventPreProcessing(){}

int NpolEventPreProcessing::GetAVNumber(const std::string &volName) {
  if(volName.substr(0,3) == "av_") {
	int underscorePos = volName.find_first_of("_",3);
	return atoi(volName.substr(3,underscorePos-3).c_str());
  } else{
	return -1;
  }
}

int NpolEventPreProcessing::GetImprNumber(const std::string &volName) {
  if(volName.substr(0,3) == "av_") {
    int underscorePos = volName.find_first_of("_",1+
											  volName.find_first_of("_",3));
    return atoi(volName.substr(underscorePos+1,1).c_str());
  } else
    return -1;
}

int NpolEventPreProcessing::GetPlacementNumber(const std::string &volName) {
  if(volName.substr(0,3) == "av_") {
    int underscorePos = volName.find_first_of("_",1+
      volName.find_first_of("_",1+
      volName.find_first_of("_",1+
      volName.find_first_of("_",1+
      volName.find_first_of("_",3)))));
    return atoi(volName.substr(underscorePos+1,std::string::npos).c_str());
  } else
    return -1;
}

void NpolEventPreProcessing::TaggerHitPosition(double hPos[],double lPos[], int detNums[]){
  double CuAnalyzerPos = 478.; // (cm) z-axis offset
  double TrackerOffset = 0;
  
  do{
	hPos[0] = lPos[0] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[0]) > 75.0);
  do{
	hPos[1] = lPos[1] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[1]) > 200.0);

  switch(detNums[1]){
  case 1:
	TrackerOffset = CuAnalyzerPos + -15.983; // (cm)
  case 2:
	TrackerOffset = CuAnalyzerPos + -9.208; // (cm)
  case 3:
	TrackerOffset = CuAnalyzerPos + 9.207; // (cm)
  case 4:
	TrackerOffset = CuAnalyzerPos + 16.547; // (cm)
  case 5:
	TrackerOffset = CuAnalyzerPos + 29.303; // (cm)
  case 6:
	TrackerOffset = CuAnalyzerPos + 42.073; // (cm)
  default:;
  }
  
  hPos[2] = lPos[2] + TrackerOffset; // keep z position the same (thin)
  RotateNpolToG4(hPos, NpolAng);  // Rotate into position of G4 System
}

void  NpolEventPreProcessing::AnalyzerHitPosition(double hPos[],double lPos[], int detNums[]){
  double AnalyzerOffset = 541.82; // (cm) z-axis offset
  
  do{
	hPos[2] = lPos[2] + rand->Gaus(0.0,1.5); //1.5cm FWHM sigma
  } while (TMath::Abs(hPos[2]) > 12.5);
  hPos[2] = hPos[2] + AnalyzerOffset;
  
  hPos[0] = rand->Uniform(0,4)-2; // (cm) randomizer
  hPos[1] = rand->Uniform(0,4)-2; // (cm)

  int rowSize = 4;
  hPos[0] = hPos[0] + (6.0 - 4.0 * float(detNums[2]%rowSize));
  hPos[1] = hPos[1] + (14.0 - 4.0 * float(detNums[2]/rowSize));

  RotateNpolToG4(hPos, NpolAng);  // Rotate into position of G4 System
   
  return;
}

//void NpolEventPreProcessing::HodoscopeHitPosition(double hPos[],double lPos[], int detNums[]){
void NpolEventPreProcessing::EarrayHitPosition(double hPos[],double lPos[], int detNums[]){
  double detYsize = 8.6; // (cm)
  double xOffset = 120.0; // (cm)
  double yOffset = detYsize *12 - 0.5*detYsize; // (cm)
  double zOffset = 495.0; // (cm) z-offset in G4 world
  
  do {
	hPos[2] = lPos[2] + rand->Gaus(0.0, 2.0); // z-axis random spread gaussian
  } while (TMath::Abs(hPos[2]) > 50.0);
  hPos[2] = hPos[2] + zOffset;
  
  hPos[0] = rand->Uniform(0,1)-0.5; // (cm)
  hPos[1] = rand->Uniform(0,8.6)-4.3; // (cm)

  hPos[1] = hPos[1] + (yOffset - detNums[2]*detYsize);

  if(detNums[1] == 1) hPos[0] = hPos[0] + xOffset;
  if(detNums[1] == 2) hPos[0] = hPos[0] - xOffset;
  
  RotateNpolToG4(hPos, NpolAng);
 
  return;
}
void NpolEventPreProcessing::DeltaEarrayHitPosition(double hPos[],double lPos[], int detNums[]){
  double CHAnalyzerPos = 541.82, HodoscopeZ = 50.0, AnalyzerZ = 25.0;  // (cm)
  double zOffset = CHAnalyzerPos + 0.5*(HodoscopeZ + AnalyzerZ);
  double xLeftHodoscope = 127.15, xRightHodoscope = 128.35;  // (cm)
  double xOffset = 0;

  do{
	hPos[1] = lPos[1] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[1]) > 200.0);
  do{
	hPos[2] = lPos[2] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[2]) > 75.0);
  hPos[2] += zOffset;

  switch(detNums[1]){
  case 1: xOffset = xLeftHodoscope - 9.63;
  case 2: xOffset = xLeftHodoscope - 16.96;
  case 3: xOffset = -(xRightHodoscope - 10.98);
  case 4: xOffset = -(xRightHodoscope - 18.31);
  default: ;
  }
  
  hPos[0] = lPos[0] + xOffset;
  RotateNpolToG4(hPos, NpolAng);
}


/*void NpolEventPreProcessing::DeltaEarrayHitPosition(double hPos[],double lPos[], int detNums[]){
  
  double VertOffSet = 42.0;
  
  do {
	hPos[0] = lPos[0] + rand->Gaus(0.0, 2.0);
  } while (TMath::Abs(hPos[0]) > 80.0);

  hPos[1] = rand->Uniform(0,1)-0.5; // (cm)
  hPos[2] = rand->Uniform(0,10)-5; // (cm)
  RotateNpolToG4(hPos, NpolAng);

  if(detNums[0] == 3) hPos[1] = hPos[1] + VertOffSet;
  if(detNums[0] == 4) hPos[1] = hPos[1] + VertOffSet + 10.0;
  if(detNums[0] == 7) hPos[1] = hPos[1] + (-VertOffSet);
  if(detNums[0] == 8) hPos[1] = hPos[1] + (-(VertOffSet + 10.0));
  
  if((detNums[0] == 3) || (detNums[0] == 7)) hPos[2] = hPos[2] + 700. + (13. - (detNums[2] + 1)) * 10.;
  if((detNums[0] == 4) || (detNums[0] == 8)) hPos[2] = hPos[2] + 830. + (14. - (detNums[2] + 1)) * 10.;
  RotateNpolToG4(hPos, NpolAng);
  return;
  }*/

/*void NpolEventPreProcessing::EarrayHitPosition(double hPos[],double lPos[], int detNums[]){
  double NDetStandardLength = 100.0;  // (cm)
  double NDetThickness = 10.0; // (cm)
  double EarrayRotAngle = 45.0 *TMath::Pi()/180.; // Erray rotation angle
  double VertOffSet = (NDetStandardLength + 60.0)/2 * sin(EarrayRotAngle) + 40.0; // 40*cm offset from geometry
  double VertOffSet2 = VertOffSet + 10.0;
  double HorOffSet = (NDetStandardLength + 60.0)/2 * cos(EarrayRotAngle) + NDetThickness/2 * sin(EarrayRotAngle);
  
  //double VertOffSet = 90.0;
  //double HorOffSet = 60.2;
  
  do {
	hPos[0] = lPos[0] + rand->Gaus(0.0, 2.0);
  } while (TMath::Abs(hPos[0]) > 80.0);
  hPos[1] = rand->Uniform(0,10)-5; // (cm)
  hPos[2] = rand->Uniform(0,10)-5; // (cm)
  RotateDetToNpol(hPos,detNums);
 
  if(((detNums[0] == 1) || (detNums[0] == 2)) && (detNums[1] == 1)) hPos[0] = hPos[0] + HorOffSet;
  if(((detNums[0] == 5) || (detNums[0] == 6)) && (detNums[1] == 1)) hPos[0] = hPos[0] + HorOffSet;
  if(((detNums[0] == 1) || (detNums[0] == 2)) && (detNums[1] == 2)) hPos[0] = hPos[0] - HorOffSet;
  if(((detNums[0] == 5) || (detNums[0] == 6)) && (detNums[1] == 2)) hPos[0] = hPos[0] - HorOffSet;
  
  RotateNpolToG4(hPos, NpolAng);
  
  if(detNums[0] == 1) hPos[1] = hPos[1] + VertOffSet;
  if(detNums[0] == 2) hPos[1] = hPos[1] + (VertOffSet2);
  if(detNums[0] == 5) hPos[1] = hPos[1] - VertOffSet;
  if(detNums[0] == 6) hPos[1] = hPos[1] - (VertOffSet2);

  if((detNums[0] == 1) || (detNums[0] == 5)) hPos[2] = hPos[2] + 700. + (13. - (detNums[2] + 1)) * 10.;
  if((detNums[0] == 2) || (detNums[0] == 6)) hPos[2] = hPos[2] + 830. + (14. - (detNums[2] + 1)) * 10.;
  
  RotateNpolToG4(hPos, NpolAng);
 
  return;
  }*/

void  NpolEventPreProcessing::RotateNpolToG4(double hPos[], double RotAng){
  // Rotate the hit point by RotAng (generally NpolAng) so hits in NPOL volumes are in the sim's G4 coordinates
  double tempPos[3] = { hPos[0], hPos[1], hPos[2] };
  double RightAng = TMath::Pi()/2;
  
  hPos[0] = TMath::Cos(RotAng) * tempPos[0] + TMath::Cos(RightAng) * tempPos[1] + TMath::Cos(RightAng+RotAng) * tempPos[2];
  hPos[1] = TMath::Cos(RightAng) * tempPos[0] + TMath::Cos(0) * tempPos[1] + TMath::Cos(RightAng) * tempPos[2];
  hPos[2] = TMath::Cos(RightAng-RotAng) * tempPos[0] + TMath::Cos(RightAng) * tempPos[1] + TMath::Cos(RotAng) * tempPos[2];

  return;
}

void  NpolEventPreProcessing::RotateDetToNpol(double hPos[], int detNums[], double detAngle){
  // A GIVEN degree rotation to rotation array hit positions into the NPOL coordinates
  double tempPos[3] = { hPos[0], hPos[1], hPos[2] };
  double RightAng = TMath::Pi()/2;
 
  if(((detNums[0] == 1) || (detNums[0] == 2)) && (detNums[1] == 2)) detAngle = -1*detAngle;
  if(((detNums[0] == 5) || (detNums[0] == 6)) && (detNums[1] == 1)) detAngle = -1*detAngle;
   
  hPos[0] = TMath::Cos(detAngle) * tempPos[0] + TMath::Cos(RightAng - detAngle) * tempPos[1] + TMath::Cos(RightAng) * tempPos[2];
  hPos[1] = TMath::Cos(RightAng + detAngle) * tempPos[0] + TMath::Cos(detAngle) * tempPos[1] + TMath::Cos(RightAng) * tempPos[2];
  hPos[2] = TMath::Cos(RightAng) * tempPos[0] + TMath::Cos(RightAng) * tempPos[1] + TMath::Cos(0) * tempPos[2];
  
  return;
}

/*void  NpolEventPreProcessing::RotateDetToNpol(double hPos[], int detNums[]){
  // A 45 degree rotation to move E-array hit positions into the NPOL coordinates
  double tempPos[3] = { hPos[0], hPos[1], hPos[2] };
  double DetAng = TMath::Pi()/4;
  double RightAng = TMath::Pi()/2;
 
  if(((detNums[0] == 1) || (detNums[0] == 2)) && (detNums[1] == 2)) DetAng = -1*DetAng;
  if(((detNums[0] == 5) || (detNums[0] == 6)) && (detNums[1] == 1)) DetAng = -1*DetAng;
   
  hPos[0] = TMath::Cos(DetAng) * tempPos[0] + TMath::Cos(RightAng - DetAng) * tempPos[1] + TMath::Cos(RightAng) * tempPos[2];
  hPos[1] = TMath::Cos(RightAng + DetAng) * tempPos[0] + TMath::Cos(DetAng) * tempPos[1] + TMath::Cos(RightAng) * tempPos[2];
  hPos[2] = TMath::Cos(RightAng) * tempPos[0] + TMath::Cos(RightAng) * tempPos[1] + TMath::Cos(0) * tempPos[2];
  
  return;
  }*/


void NpolEventPreProcessing::RotateG4ToRoot(double curPos[], double newPos[], double rotMat[][3]){ //double detAngle){
  // interchange Z and Y axis and flip sign on X to go from G4 to ROOT coordinate systems.
  //double rotMat[3][3] = { {-1.,0.,0.}, {0.,0.,1.}, {0.,1.,0.} };
  
  newPos[0] = curPos[0]*rotMat[0][0] + curPos[1]*rotMat[0][1] + curPos[2]*rotMat[0][2];
  newPos[1] = curPos[0]*rotMat[1][0] + curPos[1]*rotMat[1][1] + curPos[2]*rotMat[1][2];
  newPos[2] = curPos[0]*rotMat[2][0] + curPos[1]*rotMat[2][1] + curPos[2]*rotMat[2][2];
}
