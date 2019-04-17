
#include "NpolEventProcessing.hh"
#include "NpolEventPreProcessing.hh"
#include "NpolPhysicsVariables.hh"

#define EDEP_THRESHOLD 0.250  /*MeV*/
/*MeV This threshold is a per detector low value in SOI selection */
double NpolEventProcessing::LOW_THRESHOLD = 0.040;
/* number of analyzer layers; not general; only good for 4 and 6 layers */
int NpolEventProcessing::LAYER_NUM = 4; 


NpolEventProcessing *NpolEventProcessing::Process = NULL;

NpolEventProcessing *NpolEventProcessing::GetInstance() {
	if(Process == NULL)
		Process = new NpolEventProcessing();

	return Process;
}


NpolEventProcessing::NpolEventProcessing(){}

NpolEventProcessing::~NpolEventProcessing(){}

void NpolEventProcessing::fillDetectorEventMap(std::map<std::string,NpolDetectorEvent *> &eventMap, const NpolStep *aStep){

  // **** This method fills the detector event map which contains objects of
  // NpolDetectorEvent class keyed to the detector name (string).
  
  NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();
 
  int AVNum = PProcess->GetAVNumber(aStep->volume);
  int imprintNum = PProcess->GetImprNumber(aStep->volume);
  int PVNum = PProcess->GetPlacementNumber(aStep->volume);
    
  if(eventMap.find(aStep->volume) == eventMap.end())
	eventMap[aStep->volume] = new NpolDetectorEvent();
  
  (eventMap[aStep->volume])->totEnergyDep += aStep->eDep;
  
  if(!((eventMap[aStep->volume])->thresholdExceeded) &&
	 (eventMap[aStep->volume])->totEnergyDep >= EDEP_THRESHOLD) {
	(eventMap[aStep->volume])->thresholdExceeded = true;
	(eventMap[aStep->volume])->lPosX = aStep->lPosX;
	(eventMap[aStep->volume])->lPosY = aStep->lPosY;
	(eventMap[aStep->volume])->lPosZ = aStep->lPosZ;
	(eventMap[aStep->volume])->gPosX = aStep->gPosX;
	(eventMap[aStep->volume])->gPosY = aStep->gPosY;
	(eventMap[aStep->volume])->gPosZ = aStep->gPosZ;
	(eventMap[aStep->volume])->time = (aStep->time) + randN->Gaus(0.0,0.200);
	
	//****** Compute the hit position in the volume and save. This is done to "simulate" ****** // 
	// what we could see in a real scintillation detector based on detector resolutions.
	double hitPos[3] = { 0.0, 0.0, 0.0 };
	double lPos[3] = { aStep->lPosX, aStep->lPosY, aStep->lPosZ };
	int detNums[3] = { AVNum, imprintNum, PVNum };
	
	if((AVNum == 2)){
	  PProcess->AnalyzerHitPosition(hitPos, lPos, detNums);
	} else if((AVNum == 3)){
	  PProcess->INFNFrontTrackerHitPosition(hitPos, lPos, detNums);
	} else if((AVNum == 4)){
	  PProcess->UVaFrontTrackerHitPosition(hitPos, lPos, detNums);
	} else if((AVNum == 5)){
	  PProcess->UVaRecoilTrackerHitPosition(hitPos, lPos, detNums);
	} else if((AVNum == 6) || (AVNum == 7)){
	  PProcess->HodoscopeHitPosition(hitPos, lPos, detNums);
	}
	
	(eventMap[aStep->volume])->hPosX = hitPos[0]; 
	(eventMap[aStep->volume])->hPosY = hitPos[1]; 
	(eventMap[aStep->volume])->hPosZ = hitPos[2];

  }
}

void NpolEventProcessing::fillVertexMap(std::map<int,NpolVertex *> &theVertexMap, const std::vector<NpolVertex *> *vertVector, int DesiredPID, std::string eventVolume, std::string eventProcess, double eventTime){
  
  // **** This method fills a map that is keyed with track IDs to vertex
  // information from NpolVertex vector.  The only tracks saved in the map
  // are the original neutron (PID = 0, TID = 1) and any particle satisfying
  // DesiredPID == PID and any TID.  (PID = DesiredPID, TID = xx).
  //  Added in requirements for same volume, same event type, and same event time
  
  std::vector<NpolVertex *>::const_iterator v_it;
  for(v_it = vertVector->begin(); v_it != vertVector->end(); v_it++){
	NpolVertex *aVertex = *v_it;
	if(aVertex == NULL) continue;
	NpolVertex *copyVertex = new NpolVertex(*aVertex);
	
	int PID = aVertex->parentId;
	int TID = aVertex->trackId;
	std::string process = aVertex->process;
	std::string volName = aVertex->volume;
	double time = aVertex->time;
		
	// Well, I can't figure out how to copy data at one pointer
	// to data at another pointer (in a map) and then later delete
	// the pointers in the map without it killing the original pointers
	// So this resulted :(((
	if((PID == DesiredPID && volName == eventVolume && time == eventTime /*&& process == eventProcess*/) || PID == 0){
	  if(theVertexMap.find(TID) == theVertexMap.end()){
		theVertexMap[TID] = new NpolVertex();
		theVertexMap[TID] = copyVertex;
	  } 
	}	
  }
}

void NpolEventProcessing::printVertexMap(std::map<int,NpolVertex *> &theVertexMap, int eventID){
  
  NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();
  NpolEventProcessing *Process = NpolEventProcessing::GetInstance();
  NpolPhysicsVariables *PhysVars = NpolPhysicsVariables::GetInstance();
  
  std::cout << "*********Starting Vertex Map Dump**********" << std::endl;
  std::map<int,NpolVertex *>::iterator mapIt;
  double totalEnergy = 0.0;
  for(mapIt = theVertexMap.begin(); mapIt != theVertexMap.end(); mapIt++){
	int PID = mapIt->second->parentId;
	int TID = mapIt->second->trackId;
	int pType = mapIt->second->particleId;
	std::string partName = mapIt->second->particle;
	std::string volName = mapIt->second->volume;
	std::string process = mapIt->second->process;
	double Time = mapIt->second->time;
	double Energy = mapIt->second->energy;
	double momentum = PhysVars->computeMomentum(mapIt->second->momX, mapIt->second->momY, mapIt->second->momZ);
	if(mapIt->first != 1) totalEnergy += Energy;
	std::cout << "      Event #: " << eventID << " PID = "<< PID
			  << " TID = " << TID
			  << "   Particle " << pType << " " << partName << " AV #: "
			  << PProcess->GetAVNumber(volName)
			  << " Impr #: " << PProcess->GetImprNumber(volName) << " PV #: "
			  << PProcess->GetPlacementNumber(volName)
			  << " SOI: " << Process->sectionNumber(volName) << " Time = "
			  << Time << " Energy: " << Energy << " Process: " << process
			  << " Momentum: " << momentum << std::endl;
  }
  std::cout << " Total Energy of all particles = " << totalEnergy
			<< " MeV" << std::endl;
  std::cout << "*********Ending Vertex Map Dump**********" << std::endl;
  
}

// The NPOL polarimeter is divided into 4 or 6 sections.  This function takes
// a volume name and returns the section number that the detector belongs to,
// or -1 if the volume does not belong in the polarimeter.
int NpolEventProcessing::sectionNumber(const std::string &volName) {
  NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();
  
  int avNum, imprNum, pvNum;
  if((avNum = PProcess->GetAVNumber(volName)) == -1) return -1;

  // there is only one "secton" in RP-GEN
  return 0;

}

// This method returns the NPOL detector "type" which is defined as the following:
// analyzer = detector in 1 of the 4 vertically stacked arrays making up the analyzer sections
// tagger = detector in the arrays in front of each analyzer section
// topEArray = any detector in the top E-array detectors.  There are 2 imprints per section with 2 sections (4 total)
// botEArray = any detector in the bottom E-array detectors.  There are 2 imprints per section with 2 sections (4 total)
// topdEArray = any detector in the top dE-array detectors.  There are 2 imprints per section with 2 sections (4 total)
// botdEArray = any detector in the bottom dE-array detectors.  There are 2 imprints per section with 2 sections (4 total)
PolarimeterDetector NpolEventProcessing::detectorType(const std::string &volName) {
  NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();

  int avNum = PProcess->GetAVNumber(volName);
  int imprNum = PProcess->GetImprNumber(volName);
  switch(avNum) {
  case 1: return hCal;
  case 2: return analyzer;
  case 3: return infnGEM;
  case 4: return uvaFront;
  case 5: {
	if((imprNum == 1) || (imprNum == 2)) return topdEArray; //UVa wing GEM left
	if((imprNum == 3) || (imprNum == 4)) return botdEArray; //UVa wing GEM right
  }
  case 6: return topEArray; //Left hodoscope
  case 7: return botEArray; //Right hodoscope
  default: return unknown;
  }
}

// Return the frontmost polarimeter section that passes requirements 1 and 2.
// If -1 is returned, then no section passed requirements 1 and 2.
int NpolEventProcessing::getSectionOfInterest(const std::map<std::string,NpolDetectorEvent *> *detEvents) {
  NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();
  int sectionOfInterest = -1;
  int totalDetHit = 0; int AVNum;
  for(int section = (LAYER_NUM - 1); section >= 0; section--) {
    std::map<std::string,NpolDetectorEvent *>::const_iterator it;
	bool analyzerFlag = false;
	bool topdEArrayFlag = false;
	bool botdEArrayFlag = false;
	bool taggerFlag = false;
	for(it = detEvents->begin(); it != detEvents->end(); it++) {
	  if(sectionNumber(it->first) == section) {
		AVNum =  PProcess->GetAVNumber(it->first);
		if((AVNum >= 2) && (AVNum <=7)){
		  if(it->second->totEnergyDep >= LOW_THRESHOLD) totalDetHit++; 
		}
		switch(detectorType(it->first)) {
		case analyzer: analyzerFlag |= it->second->thresholdExceeded == true;
		  break;
		case tagger: taggerFlag |= it->second->thresholdExceeded == true; 
		  break;
		case topdEArray: topdEArrayFlag |= it->second->thresholdExceeded == true; 
		  break;
		case botdEArray: botdEArrayFlag |= it->second->thresholdExceeded == true; 
		  break;
		default:
		  break;
		}
	  }
	}
	if( ((analyzerFlag) /*&& (!taggerFlag)*/) && ((topdEArrayFlag) || (botdEArrayFlag)) ) sectionOfInterest = section;
  }

  // With the preliminaary selection of SOI, now check the E-arrays in the SOI and SOI+1 for at least 1 hit
  if(sectionOfInterest != -1){
	int firstSection = sectionOfInterest;
	int secondSection;
	bool topEArrayFlag = false;
	bool botEArrayFlag = false;
	if(sectionOfInterest == (LAYER_NUM - 1)){
	  secondSection = LAYER_NUM - 1;
	} else {
	  secondSection = sectionOfInterest + 1;
	}
	for(int section = firstSection; section <= secondSection; section++) {
	  std::map<std::string,NpolDetectorEvent *>::const_iterator E_it;
	  for(E_it = detEvents->begin(); E_it != detEvents->end(); E_it++) {
		if(sectionNumber(E_it->first) == section) {
		  switch(detectorType(E_it->first)) {
		  case topEArray: topEArrayFlag |= E_it->second->thresholdExceeded == true; 
			break;
		  case botEArray: botEArrayFlag |= E_it->second->thresholdExceeded == true; 
			break;
		  default:
			break;
		  }
		}
	  }
	}
	if(!((topEArrayFlag) || (botEArrayFlag))) sectionOfInterest = -1;
  }

  // Reject events with more than 40 detectors with hits above 0.040 MeV; Otherwise return the ID'd section of interest
  /*if((sectionOfInterest != -1) && (totalDetHit >= 40)) { 
	std::cout << "Event #  Rejected! Total number of detectors with 40 keV or greater: " << totalDetHit << std::endl;
	  sectionOfInterest = -1;
	  }*/
  
  return sectionOfInterest;
}

PolarimeterDetector NpolEventProcessing::getEArrayOfInterest(std::map<std::string,NpolDetectorEvent *> *detEvents, std::map<PolarimeterDetector, double> *eDepArrayTotal, int SOI) {
  double topETotal = (*eDepArrayTotal)[topEArray];
  double botETotal = (*eDepArrayTotal)[botEArray];
  double topdETotal = (*eDepArrayTotal)[topdEArray];
  double botdETotal = (*eDepArrayTotal)[botdEArray];

  PolarimeterDetector EOI = unknown;
  double totalTop = topETotal + topdETotal;
  double totalBot = botETotal + botdETotal;
  double topR = totalTop/totalBot;
  double botR = totalBot/totalTop;
  bool topEFlag = false;
  bool botEFlag = false;
  bool topdEFlag = false;
  bool botdEFlag = false;
  
  if((topETotal >= 5) && (topETotal <= 150)) topEFlag = true;  // energy units of MeV
  if((botETotal >= 5) && (botETotal <= 150)) botEFlag = true;  // energy units of MeV
  if((topdETotal >= 1) && (topdETotal <= 25)) topdEFlag = true;  // energy units of MeV
  if((botdETotal >= 1) && (botdETotal <= 25)) botdEFlag = true;  // energy units of MeV
  
  if(topEFlag && /*(topR >= 20) &&*/ !(botEFlag)){
    if(topdEFlag && !(botdEFlag)) EOI = topEArray;
  } else if(!(topEFlag) && botEFlag /*&& (botR >= 20)*/){
    if(!(topdEFlag) && botdEFlag) EOI = botEArray;
  } else {
    EOI = unknown;
	//std::cout << "Failed to detect top or bottom E-array" << std::endl;
  }

  return EOI;
}

void NpolEventProcessing::getEDepArrayTotal(std::map<std::string,NpolDetectorEvent *> *detEvents, std::map<PolarimeterDetector, double> *eDepArrayTotal, int SOI){
  NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();
  
  std::map<std::string,NpolDetectorEvent *>::iterator e_it;
  for(e_it = detEvents->begin(); e_it != detEvents->end(); e_it++) {	
	int section = sectionNumber(e_it->first);
	PolarimeterDetector detector = detectorType(e_it->first);
	if(section == SOI) {
	  if((detector == analyzer) || (detector == tagger) || (detector == topdEArray)
		 || (detector == botdEArray) || (detector == topEArray) || (detector == botEArray)){
		(*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
	  }
	}


	
	/*if(section == (SOI + 1)){
	  int AVNum = PProcess->GetAVNumber(e_it->first);
	  int PVNum = PProcess->GetPlacementNumber(e_it->first);
	  if(section == 1){
		if((detector == topdEArray) || (detector == botdEArray)){
		  if((AVNum == 3) || (AVNum == 7)){
			if((PVNum == 5)){
			  (*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
			}
		  }
		} else if((detector == topEArray) || (detector == botEArray)){
		  if((AVNum == 1) || (AVNum == 5)){
			if((PVNum >= 3) && (PVNum <= 5)){
			  (*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
			} 
		  }
		}
	  } else if(section == 2){
		if((detector == topdEArray) || (detector == botdEArray)){
		  if((AVNum == 4) || (AVNum == 8)){
			if((PVNum == 13)){
			  (*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
			}
		  }
		} else if((detector == topEArray) || (detector == botEArray)){
		  if((AVNum == 2) || (AVNum == 6)){
			if((PVNum >= 11) && (PVNum <= 13)){
			  (*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
			} 
		  }
		}
	  } else if(section == 3){
		if((detector == topdEArray) || (detector == botdEArray)){
		  if((AVNum == 4) || (AVNum == 8)){
			if((PVNum == 6)){
			  (*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
			}
		  }
		} else if((detector == topEArray) || (detector == botEArray)){
		  if((AVNum == 2) || (AVNum == 6)){
			if((PVNum >= 4) && (PVNum <= 6)){
			  (*eDepArrayTotal)[detector] += e_it->second->totEnergyDep;
			} 
		  }
		}
	  }
	  }*/
  }
  return;
}

bool NpolEventProcessing::checkEarrayHits(const std::map<std::string,NpolDetectorEvent *> *detEvents){
  std::map<std::string,NpolDetectorEvent *>::const_iterator it;
  int countTop = 0;
  int countBottom = 0;
  for(it = detEvents->begin(); it != detEvents->end(); it++) {
	if((Process->detectorType(it->first) == topEArray) && (it->second->thresholdExceeded == true)) countTop++;
	if((Process->detectorType(it->first) == botEArray) && (it->second->thresholdExceeded == true)) countBottom++;
  }
  std::cout << "Number Top E-Array Detectors 'Hit' =  " << countTop << std::endl;
  std::cout << "Number Bottom E-Array Detectors 'Hit' =  " << countBottom << std::endl;
  return true;
}

bool NpolEventProcessing::checkdEarrayHits(const std::map<std::string,NpolDetectorEvent *> *detEvents){
  std::map<std::string,NpolDetectorEvent *>::const_iterator it;
  int countTop = 0;
  int countBottom = 0;
  for(it = detEvents->begin(); it != detEvents->end(); it++) {
	if((Process->detectorType(it->first) == topdEArray) && (it->second->thresholdExceeded == true)) countTop++;
	if((Process->detectorType(it->first) == botdEArray) && (it->second->thresholdExceeded == true)) countBottom++;
  }
  std::cout << "Number Top dE-Array Detectors 'Hit' =  " << countTop << std::endl;
  std::cout << "Number Bottom dE-Array Detectors 'Hit' =  " << countBottom << std::endl;
  return true;
}

void NpolEventProcessing::fillFourVector(const NpolStep *aStep, std::pair<double,std::vector<double> > &a4Vector){

  std::vector<double> a3Vector = {0.0, 0.0,0.0};
  double aEnergy = 0.0;
  if(aStep != NULL) {
	aEnergy = aStep->energy;
	a3Vector[0] = aStep->momX;
	a3Vector[1] = aStep->momY;
	a3Vector[2] = aStep->momZ;
	a4Vector.first = aEnergy;
	a4Vector.second = a3Vector;
  } 
  a4Vector.first = aEnergy;
  a4Vector.second = a3Vector;
}

void NpolEventProcessing::fillFourVector(const NpolVertex *aVertex, std::pair<double,std::vector<double> > &a4Vector){
  std::vector<double> a3Vector = {0.0, 0.0,0.0};
  double aEnergy = 0.0;
  if(aVertex != NULL){
	aEnergy = aVertex->energy;
	a3Vector[0] = aVertex->momX;
	a3Vector[1] = aVertex->momY;
	a3Vector[2] = aVertex->momZ;
  }
  a4Vector.first = aEnergy;
  a4Vector.second = a3Vector;
}

void NpolEventProcessing::fillFourVector(std::vector<double> &particleData, std::pair<double,std::vector<double> > &a4Vector){
  
  std::vector<double> a3Vector = {0.0, 0.0,0.0};
  double aEnergy = particleData[0];
  a3Vector[0] = particleData[1];
  a3Vector[1] = particleData[2];
  a3Vector[2] = particleData[3];
  a4Vector.first = aEnergy;
  a4Vector.second = a3Vector;

}
