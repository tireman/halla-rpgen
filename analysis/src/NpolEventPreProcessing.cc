
#include "NpolEventPreProcessing.hh"

double NpolEventPreProcessing::AnalyzerZ = 25.0; // (cm)
double NpolEventPreProcessing::HodoscopeZ = 50.0; // (cm)
double NpolEventPreProcessing::NpolAng = 0.431096; /*radians; angle of NPOL rel. to beam axis*/
double NpolEventPreProcessing::CuAnalyzerPos = 478.0; // (cm) to the center of Cu Analyzer
double NpolEventPreProcessing::CHAnalyzerPos = CuAnalyzerPos + 63.82; // (cm) offset from Cu
double NpolEventPreProcessing::LeftHodoXPos = 127.15; // (cm) x-direction (left right of beam)
double NpolEventPreProcessing::RightHodoXPos = 128.35; // (cm) x-direction (left right of beam)
double NpolEventPreProcessing::HodoYPos = 0.0;  // (cm) y-direction (up-down of beam)
double NpolEventPreProcessing::HodoZPos = CHAnalyzerPos + 0.5*(HodoscopeZ + AnalyzerZ);


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

void NpolEventPreProcessing::UVaFrontTrackerHitPosition(double hPos[],double lPos[], int detNums[]){
  double TrackerOffsetZ = 0;
  
  do{
	hPos[0] = lPos[0] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[0]) > 37.5);
  do{
	hPos[1] = lPos[1] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[1]) > 100.0);

  switch(detNums[1]){
  case 1: TrackerOffsetZ = CuAnalyzerPos + -15.983; // (cm)
  case 2: TrackerOffsetZ = CuAnalyzerPos + -9.208; // (cm)
  case 3: TrackerOffsetZ = CuAnalyzerPos + 9.207; // (cm)
  case 4: TrackerOffsetZ = CuAnalyzerPos + 16.547; // (cm)
  case 5: TrackerOffsetZ = CuAnalyzerPos + 29.303; // (cm)
  case 6: TrackerOffsetZ = CuAnalyzerPos + 42.073; // (cm)
  default: ;
  }

  hPos[2] = rand->Uniform(0,0.250)-0.125; // (cm) random position in z-axis of tracker plane
  hPos[2] += TrackerOffsetZ; // keep z position the same (thin)
  RotateNpolToG4(hPos, NpolAng);  // Rotate into position of G4 System
}

void NpolEventPreProcessing::INFNFrontTrackerHitPosition(double hPos[],double lPos[], int detNums[]){
  double TrackerOffsetZ = 0;
  
  do{
	hPos[0] = lPos[0] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[0]) > 20.0);
  do{
	hPos[1] = lPos[1] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[1]) > 60.0);

  switch(detNums[1]){
  case 1: TrackerOffsetZ = CuAnalyzerPos + -39.53; // (cm)
  case 2: TrackerOffsetZ = CuAnalyzerPos + -27.24; // (cm)
  default: ;
  }

  hPos[2] = rand->Uniform(0,0.250)-0.125; // (cm) random position in z-axis of tracker plane
  hPos[2] += TrackerOffsetZ; // keep z position the same (thin)
  RotateNpolToG4(hPos, NpolAng);  // Rotate into position of G4 System
}

void  NpolEventPreProcessing::AnalyzerHitPosition(double hPos[],double lPos[], int detNums[]){
  
  do{ 
	hPos[2] = lPos[2] + rand->Gaus(0.0,1.5); //1.5cm FWHM sigma
  } while (TMath::Abs(hPos[2]) > 12.5);
  hPos[2] = hPos[2] + CHAnalyzerPos;
  
  hPos[0] = rand->Uniform(0,4)-2; // (cm) randomizer
  hPos[1] = rand->Uniform(0,4)-2; // (cm)

  int rowSize = 4;
  hPos[0] = hPos[0] + (6.0 - 4.0 * float(detNums[2]%rowSize));
  hPos[1] = hPos[1] + (14.0 - 4.0 * float(detNums[2]/rowSize));

  RotateNpolToG4(hPos, NpolAng);  // Rotate into position of G4 System
   
  return;
}

void NpolEventPreProcessing::HodoscopeHitPosition(double hPos[],double lPos[], int detNums[]){
  double detYsize = 8.6; // (cm)
  double yOffset = detYsize *12 - 0.5*detYsize; // (cm)

  // z and y dimensions the same on both hodoscopes; x is different
  do {
	hPos[2] = lPos[2] + rand->Gaus(0.0, 2.0); // z-axis random spread gaussian
  } while (TMath::Abs(hPos[2]) > 25.0);
  hPos[2] += HodoZPos;
  
  hPos[1] = rand->Uniform(0,8.6)-4.3; // (cm)
  hPos[1] += HodoYPos + (yOffset - detNums[2]*detYsize);
  
  // left hodoscope is 3 mm thick (x dimension)
  if(detNums[0] == 6) {
	hPos[0] = rand->Uniform(0,0.3)-0.15; // (cm)
	hPos[0] += LeftHodoXPos;
  }
  
  // left hodoscope is 3 cm thick (x dimension)
  if(detNums[0] == 7) {
	hPos[0] = rand->Uniform(0,3.0)-1.5; // (cm)
	hPos[0] += -RightHodoXPos;
  }
  
  RotateNpolToG4(hPos, NpolAng);
  return;
}

void NpolEventPreProcessing::UVaRecoilTrackerHitPosition(double hPos[],double lPos[], int detNums[]){
  double xOffset = 0;

  do{
	hPos[1] = lPos[1] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[1]) > 100.0);
  do{
	hPos[2] = lPos[2] + rand->Gaus(0.0,0.0070); // 70*um FWHM sigma
  } while (TMath::Abs(hPos[2]) > 37.5);
  hPos[2] += HodoZPos;

  switch(detNums[1]){
  case 1: xOffset = LeftHodoXPos - 9.63;
  case 2: xOffset = LeftHodoXPos - 16.96;
  case 3: xOffset = -(RightHodoXPos - 10.98);
  case 4: xOffset = -(RightHodoXPos - 18.31);
  default: ;
  }
  
  hPos[0] = rand->Uniform(0,0.250)-0.125; // (cm)
  hPos[0] += xOffset;
  RotateNpolToG4(hPos, NpolAng);
}


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

void NpolEventPreProcessing::RotateG4ToRoot(double curPos[], double newPos[], double rotMat[][3]){ 
  // interchange Z and Y axis and flip sign on X to go from G4 to ROOT coordinate systems.
   
  newPos[0] = curPos[0]*rotMat[0][0] + curPos[1]*rotMat[0][1] + curPos[2]*rotMat[0][2];
  newPos[1] = curPos[0]*rotMat[1][0] + curPos[1]*rotMat[1][1] + curPos[2]*rotMat[1][2];
  newPos[2] = curPos[0]*rotMat[2][0] + curPos[1]*rotMat[2][1] + curPos[2]*rotMat[2][2];
}
