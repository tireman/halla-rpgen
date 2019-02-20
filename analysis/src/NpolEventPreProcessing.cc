

#include "NpolEventPreProcessing.hh"

double NpolEventPreProcessing::NpolAng = 0.431096; /*radians; 24.7 deg angle of NPOL relative to beam axis*/

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

void  NpolEventPreProcessing::AnalyzerHitPosition(double hPos[],double lPos[], int detNums[]){
  double AnalyzerOffset = 470.; // (cm) z-axis offset
  
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

void  NpolEventPreProcessing::BackTaggerHitPosition(double hPos[],double lPos[], int detNums[]){

  double xyOffset = 75.0; // (cm)
  double zOffset = 565.0; //(cm) tagger position in G4 world

  do {
	hPos[0] = lPos[0] + rand->Gaus(0.0, 2.0); // randomize x-axis spread
  } while (TMath::Abs(hPos[0]) > 80.0);
  
  hPos[1] = rand->Uniform(0,10)-5; // (cm) y-axis random position
  hPos[2] = rand->Uniform(0,1)-0.5; // (cm) z-axis random position
  
  hPos[1] = hPos[1] + (xyOffset - 10.0 * float(detNums[2]));  // offset in vertical direction


  if(detNums[1] == 1){  // first layer is horizontal to floor

	hPos[2] = hPos[2] + zOffset;  // offset along z-axis

  } else if (detNums[1] == 2){  // second layer is 90-deg to first (vertical)
	
	hPos[2] = hPos[2] + zOffset + 2.0;  // offset along z-axis plus 2 cm as imprint 2 is farther away
	double detAngle = -1*TMath::Pi()/2; // counter clockwise rotation
	RotateDetToNpol(hPos,detNums,detAngle);
  }
  
  RotateNpolToG4(hPos, NpolAng);
  return;
}

void NpolEventPreProcessing::HodoscopeHitPosition(double hPos[],double lPos[], int detNums[]){

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

void NpolEventPreProcessing::RotateG4ToRoot(double curPos[], double newPos[], double detAngle){
  // interchange Z and Y axis and flip sign on X to go from G4 to ROOT coordinate systems.
  double rotMat[3][3] = { {-1.,0.,0.}, {0.,0.,1.}, {0.,1.,0.} };
  
  newPos[0] = curPos[0]*rotMat[0][0] + curPos[1]*rotMat[0][1] + curPos[2]*rotMat[0][2];
  newPos[1] = curPos[0]*rotMat[1][0] + curPos[1]*rotMat[1][1] + curPos[2]*rotMat[1][2];
  newPos[2] = curPos[0]*rotMat[2][0] + curPos[1]*rotMat[2][1] + curPos[2]*rotMat[2][2];
}
