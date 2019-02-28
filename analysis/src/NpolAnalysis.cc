/* Npol Analysis Script is designed to analyze the neutron flux on 
   the NPOL polarimeter being designed by the CGEN collaboration at 
   Jefferson National Laboratory. (2016) Revisions: Spring 2017 by 
   Will Tireman and Ashley Adzima (added some histograms) Revisions: 
   January-March 2018 by Will Tireman (Fixed eff. calculation, cleaned up code)
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>

#include "TVector3.h"
#include <TVectorD.h>
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TRandom3.h"

#include "NpolVertex.hh"
#include "NpolTagger.hh"
#include "NpolStep.hh"
#include "NpolStatistics.hh"
#include "NpolDetectorEvent.hh"

#include "NpolFileEnvManager.hh"
#include "NpolEventPreProcessing.hh"
#include "NpolEventProcessing.hh"
#include "NpolPhysicsVariables.hh"
#include "NpolHistoManager.hh"

using namespace std;

//********************** Definition of Variables and Classes *********************//

#define EDEP_THRESHOLD 1.0  /*MeV*/
#define angleLow 45.3       /*degrees: low angle recoil proton cut*/
#define angleHigh 81.6      /*degrees; high angle recoil proton cut*/

NpolEventProcessing *Process = NpolEventProcessing::GetInstance();
NpolEventPreProcessing *PProcess = NpolEventPreProcessing::GetInstance();
NpolPhysicsVariables *PhysVars = NpolPhysicsVariables::GetInstance();
NpolFileEnvManager *FEman = NpolFileEnvManager::GetInstance();
NpolHistoManager *HistoMan = NpolHistoManager::GetInstance();

//***************** End Definitions of Variables and Classes ***********************//

//************************ Main Program Begins Here ********************************//

int main(int argc, char *argv[]) {
  
  FEman->RetrieveENVvariables();

  TString InputFile = FEman->FormInputFile(FEman->InputDir);
  TString OutputFile = FEman->FormOutputFile(FEman->OutputDir);

  HistoMan->OpenFile(OutputFile);

  TChain *npolTree = new TChain("T");
  TChain *statsTree = new TChain("T2");

  npolTree->Add(InputFile);
  statsTree->Add(InputFile);

  std::vector<NpolStep *> *steps = NULL;
  std::vector<NpolVertex *> *verts = NULL;
  std::vector<NpolTagger *> *tagEvent = NULL;
  std::vector<NpolStatistics *> *stats = NULL;
  npolTree->SetBranchAddress("steps",&steps);
  npolTree->SetBranchAddress("tracks",&verts);
  npolTree->SetBranchAddress("NPOLTagger",&tagEvent);
  statsTree->SetBranchAddress("stats",&stats);
  
  //********************************* Define your Histograms Here *******************************

  //    1D-Histograms
  HistoMan->CreateHistograms("Analyzer","Analyzer hit counts",34,20999,21033);
  HistoMan->CreateHistograms("HCAL","HCAL hit counts",290,10999,11289); 
  HistoMan->CreateHistograms("LeftHodo","Left Hodoscope hit counts",26,60999,61025);
  HistoMan->CreateHistograms("RightHodo","Right Hodoscope hit counts",26,70999,71025);

  for(int i = 0; i < 288; i++){
	std::string hName = "HCAL_Edep_" + std::to_string(i+1);
	std::string hTitle = "HCAL Detector " + std::to_string(i+1) + " Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,200);
  }
  for(int i = 0; i < 32; i++){
	std::string hName = "Analyzer_Edep_" + std::to_string(i+1);
	std::string hTitle = "Analyzer Detector " + std::to_string(i+1) + " Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }
  for(int i = 0; i < 24; i++){
	std::string hName = "LHodo_Edep_" + std::to_string(i+1);
	std::string hTitle = "Left Hodoscope Detector " + std::to_string(i+1) + " Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
	hName = "RHodo_Edep_" + std::to_string(i+1);
	hTitle = "Right Hodoscope Detector " + std::to_string(i+1) + " Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }
  for(int i = 0; i < 2; i++){
	std::string hName = "INFN_Front_Edep_" + std::to_string(i+1);
	std::string hTitle = "INFN Front Tracker GEM Detector " + std::to_string(i+1) + "Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }
  for(int i = 0; i < 2; i++){
	std::string hName = "UVA_Front_Edep_" + std::to_string(i+1);
	std::string hTitle = "UVA Front Tracker GEM Detector " + std::to_string(i+1) + "Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }
  for(int i = 0; i < 4; i++){
	std::string hName = "UVA_Rear_Edep_" + std::to_string(i+1);
	std::string hTitle = "UVA Rear Tracker GEM Detector " + std::to_string(i+1) + "Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }
  for(int i = 0; i < 2; i++){
	std::string hName = "UVA_LeftWing_Edep_" + std::to_string(i+1);
	std::string hTitle = "UVA Left Wing Tracker GEM Detector " + std::to_string(i+1) + "Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }
  for(int i = 0; i < 2; i++){
	std::string hName = "UVA_RightWing_Edep_" + std::to_string(i+1);
	std::string hTitle = "UVA Right Wing Tracker GEM Detector " + std::to_string(i+1) + "Energy Deposited (MeV)";
	HistoMan->CreateHistograms(hName.c_str(),hTitle.c_str(),100,0,20);
  }

  //    2D-Histograms
  //HistoMan->CreateHistograms("Generic_2D_Histo","A Generic 2D Histogram",200,-350.0,350.0, 200,-50.0,50.0);
  
  //    3D-Histograms
  //HistoMan->CreateHistograms("Hit_Position_3D","3D Hit Position Histogram",200, 0.0,550.0, 200,350.0,575.0, 200,-125.0,125.0); 
  //HistoMan->CreateHistograms("Global_Position_3D","3D Global Hit Position Histogram",200,0.0,550.0, 200,350.,575., 200,-125,125); 
  //********************************* End Histogram Definitions ********************************

 
  // ****** BEGIN STATS LOOP ******
  int totalEvents = 0;   // Total number of neutrons (events) generated and recorded in the data file
  int taggedEvents = 0;  // Total number of neutrons (events) which make it through magnets and lead curtain
  int eventsPassed = 0;  // Total number of neutron scattering events which pass all Proposal 37 cuts
  int eventsFailed = 0;  // Total number of neutron scattering events which do not pass all Proposal 37 cuts
  for(int i = 0; i < statsTree->GetEntries(); i++) {
    statsTree->GetEntry(i);
    totalEvents += ((*stats)[0])->totalEvents;
  }
  double electronTime = totalEvents/(6.242e12); //6.242e12 e-/s at 1 microAmp (amount of physical time at 1uA)
  std::cout << totalEvents << " electrons thrown at setup." << std::endl;
  // ****** END STATS LOOP ****** 
  
  std::map<int,int> scintCounter;
  // ****** BEGIN EVENT LOOP ****** 
  int nentries = npolTree->GetEntries();
  TRandom *rand = new TRandom3();
  for(int i = 0; i < nentries; i++) {
  //for(int i = 0; i < 1000; i++) {
   	if(i % 1000 == 0) std::cout << "Processing event #" << i << std::endl;
	npolTree->GetEntry(i);
	
	//std::pair<double,std::vector<double> > initNeutron4Vec;
	//std::pair<double,std::vector<double> > recoilParticle4Vec;
	//std::pair<double,std::vector<double> > projNeutron4Vec;
	//std::pair<double,std::vector<double> > scattNeutron4Vec;
	//std::pair<double,std::vector<double> > scattParticle4Vec;

	std::map<std::string,NpolVertex *> vertexMap;  // Vertex Map Keyed to volume of origin
	std::map<std::string,NpolDetectorEvent *> detEvents;   // Event map (NPOL Detector Class)
	std::map<PolarimeterDetector, double> eDepArrayTotal;  // Total energy map for each array
    eDepArrayTotal[analyzer] = 0.0;
    eDepArrayTotal[tagger] = 0.0;
    eDepArrayTotal[topEArray] = 0.0;
    eDepArrayTotal[topdEArray] = 0.0;
    eDepArrayTotal[botEArray] = 0.0;
    eDepArrayTotal[botdEArray] = 0.0;

	// Cycle through the tracks vector and store information in a map
	std::vector<NpolVertex *>::iterator v_it;
    for(v_it = verts->begin(); v_it != verts->end(); v_it++){
      NpolVertex *vertex = *v_it;
      if(vertex == NULL) continue;
	  NpolVertex *copyVertex = new NpolVertex(*vertex);
	  
	  std::string volumeName = vertex->volume;
	  if(vertexMap.find(volumeName) == vertexMap.end()){
		vertexMap[volumeName] = new NpolVertex();
		vertexMap[volumeName] = copyVertex;
	  }
	}
	// End the tracks loop

	// fill detector map (as defined in NpolLib) with steps
	std::vector<NpolStep *>::iterator s_it;
	std::vector<NpolTagger *>::iterator t_it;
	bool eventFlag = false;
    for(s_it = steps->begin(); s_it != steps->end(); s_it++) {
      NpolStep *aStep = *s_it;
	  if(aStep == NULL) continue;
	  Process->fillDetectorEventMap(detEvents,aStep);
	}
	// End steps loop


	// Cycle through the detector event map and do stuff
	//double curPos[3];
	//double newPos[3];
	//double rotMat[3][3] = { {-1.,0.,0.}, {0.,0.,1.}, {0.,1.,0.} };
	std::map<std::string,NpolDetectorEvent *>::iterator det_it;
	for(det_it = detEvents.begin(); det_it != detEvents.end(); det_it++) {
	  std::string volumeName = det_it->first;
	  int AVNum = PProcess->GetAVNumber(volumeName);
	  int ImprNum = PProcess->GetImprNumber(volumeName);
	  int PVNum = PProcess->GetPlacementNumber(volumeName);
	  double eDep = 0.0;
	  std::string hName = "";
	  if(AVNum != -1 && det_it->second->thresholdExceeded){
		int detNum = AVNum*10000 + ImprNum *1000 + PVNum;
		if(scintCounter.find(detNum) == scintCounter.end()) scintCounter[detNum] = 0;
		scintCounter[detNum]++;
	  }
	 
	  eDep = det_it->second->totEnergyDep;
	  if(eDep == 0) continue;
	  if (AVNum == 1){
		hName = "HCAL_Edep_" + std::to_string(PVNum + 1);
		HistoMan->FillHistograms(hName,eDep);
	  } else if(AVNum == 2){
		hName = "Analyzer_Edep_" + std::to_string(PVNum + 1);
		HistoMan->FillHistograms(hName,eDep);
	  } else if (AVNum == 3){
		if(ImprNum == 1) hName = "INFN_Front_Edep_1";
		if(ImprNum == 2) hName = "INFN_Front_Edep_2";
		HistoMan->FillHistograms(hName,eDep);
	  } else if (AVNum == 4){
		if(ImprNum == 1) hName = "UVA_Front_Edep_1";
		if(ImprNum == 2) hName = "UVA_Front_Edep_2";
		if(ImprNum == 3) hName = "UVA_Rear_Edep_1";
		if(ImprNum == 4) hName = "UVA_Rear_Edep_2";
		if(ImprNum == 5) hName = "UVA_Rear_Edep_3";
		if(ImprNum == 6) hName = "UVA_Rear_Edep_4";
		HistoMan->FillHistograms(hName,eDep);
	  } else if (AVNum == 5){
		if(ImprNum == 2) hName = "UVA_LeftWing_Edep_1";
		if(ImprNum == 1) hName = "UVA_LeftWing_Edep_2";
		if(ImprNum == 4) hName = "UVA_RightWing_Edep_1";
		if(ImprNum == 3) hName = "UVA_RightWing_Edep_2";
		HistoMan->FillHistograms(hName,eDep);
	  } else if (AVNum == 6){
		hName = "LHodo_Edep_" + std::to_string(PVNum + 1);
		HistoMan->FillHistograms(hName,eDep);
	  } else if (AVNum == 7){
		hName = "RHodo_Edep_" + std::to_string(PVNum + 1);
		HistoMan->FillHistograms(hName,eDep);
	  }
	}
	 	
    // Clear out the maps for the next event
	
	eDepArrayTotal.clear();
	std::map<std::string,NpolDetectorEvent *>::iterator e_it;
    for(e_it = detEvents.begin(); e_it != detEvents.end(); e_it++) delete e_it->second;
    detEvents.clear();
	
	std::map<std::string,NpolVertex *>::iterator v_it2;
	for(v_it2 = vertexMap.begin(); v_it2 != vertexMap.end(); v_it2++) delete v_it2->second;
	vertexMap.clear();
	
	//initNeutron4Vec.second.clear();
	//recoilParticle4Vec.second.clear();
	//projNeutron4Vec.second.clear();
	//scattNeutron4Vec.second.clear();
	//scattParticle4Vec.second.clear();
	
  }	// END EVENT LOOP
  
  std::map<int,int>::iterator scint_it;
  for(scint_it = scintCounter.begin(); scint_it != scintCounter.end(); scint_it++){
	if(scint_it->first >= 60099 && scint_it->first <= 61025) {
	  for(int j = 0; j < scint_it->second; j++){
		HistoMan->FillHistograms("LeftHodo",scint_it->first);
	  }
	} else if(scint_it->first >= 70099 && scint_it->first <= 71025) {
	  for(int j = 0; j < scint_it->second; j++){
		HistoMan->FillHistograms("RightHodo",scint_it->first);
	  }
	} else if(scint_it->first >= 20099 && scint_it->first <= 21033) {
	  for(int j = 0; j < scint_it->second; j++){
		HistoMan->FillHistograms("Analyzer",scint_it->first);
	  }
	} else if(scint_it->first >= 11000 && scint_it->first <= 11289) {
	  for(int j = 0; j < scint_it->second; j++){
		HistoMan->FillHistograms("HCAL",scint_it->first);
	  }
	}
  }

  /*double counts = 0;
  std::map<std::string,TH1F *>::iterator it_1D;
  for(it_1D == histoMap.begin(); it_1D < histoMap.end(); it_1D++){
	counts = it_1D->second->Integral();
	std::cout << "Counts = " << counts << std::endl;
	}*/
  
  // Print some stats to the terminal window
  std::cout << taggedEvents << " of the initial " << totalEvents << " neutrons have crossed the Npol Tagger."
			<< std::endl;
  std::cout << eventsPassed << " events passed requirements.  "
			<< (taggedEvents - eventsPassed) << " failed." << std::endl;
  std::cout << (((double)eventsPassed)*100)/(double)taggedEvents << " % of the " << taggedEvents
			<< " neutrons passed cuts." << std::endl;
  
  
  // **** Fill the statistics vector ****** //
  TVectorD runStatistics(6);
  runStatistics[0] = static_cast<double>(totalEvents);
  runStatistics[1] = static_cast<double>(taggedEvents);
  runStatistics[2] = static_cast<double>(eventsPassed);
  runStatistics[3] = static_cast<double>(eventsFailed);
  runStatistics[4] = static_cast<double>(-1.0);         
  runStatistics[5] = static_cast<double>(-1.0);
  
  // ****** Write out all the results to the root file and close the file(s) ******** //
  runStatistics.Write();
  HistoMan->WriteHistograms();
  HistoMan->CloseFile();
  return 0;
}

// ************** End Main Program here *************** //

// ************* Methods for Main Program are below  this line ************** //

// No methods here ... they have all been moved to separate classes!



  //npolTree->SetCacheSize(50000000);  // This increases the amount of data loaded 
  //statsTree->SetCacheSize(50000000); // per call for more data to chew on. 
  
	//HistoMan->FillHistograms("Generic_2D_Histo",newPos[0],newPos[2]);
		
		/*curPos[0] = det_it->second->gPosX;
		curPos[1] = det_it->second->gPosY;
		curPos[2] = det_it->second->gPosZ;
		PProcess->RotateG4ToRoot(curPos,newPos,0);
		HistoMan->FillHistograms("Global_Position_3D",newPos[0], newPos[1], newPos[2]);

		curPos[0] = det_it->second->hPosX;
		curPos[1] = det_it->second->hPosY;
		curPos[2] = det_it->second->hPosZ;
		PProcess->RotateG4ToRoot(curPos,newPos,0);
		HistoMan->FillHistograms("Hit_Position_3D",newPos[0], newPos[1], newPos[2]);*/
