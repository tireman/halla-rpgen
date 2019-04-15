
#include "NpolHistoManager.hh"

NpolHistoManager *NpolHistoManager::HistoMan = NULL;

NpolHistoManager *NpolHistoManager::GetInstance() {
  if(HistoMan == NULL) HistoMan = new NpolHistoManager();
  return HistoMan;
}

NpolHistoManager::NpolHistoManager(){}

NpolHistoManager::~NpolHistoManager(){
  DeleteHistograms();
  ClearHistograms();
  delete outFile2;
}

void NpolHistoManager::CreateHistograms(std::string hID, std::string hTitle,int xBins, double xMin, double xMax){
  if(xBins == 0){
	std::cout << "X Bins value must be a non-zero integer" << std::endl;
	return;
  }

  std::string histogramRef = "h_" + hID;
  histoMap[histogramRef] = new TH1F(hID.c_str(),hTitle.c_str(),xBins,xMin,xMax);
}



void NpolHistoManager::CreateHistograms(std::string hID, std::string hTitle,int xBins, double xMin, double xMax, int yBins, double yMin, double yMax){
  if(xBins == 0){
	std::cout << "X Bins value must be a non-zero integer" << std::endl;
	return;
  } else if(yBins == 0){
	std::cout << "Y Bins value must be a non-zero integer" << std::endl;
	return;
  }

  std::string histogramRef = "h_" + hID;
  histoMap2D[histogramRef] = new TH2F(hID.c_str(),hTitle.c_str(),xBins,xMin,xMax,yBins,yMin,yMax);
}

void NpolHistoManager::CreateHistograms(std::string hID, std::string hTitle,int xBins, double xMin, double xMax, int yBins, double yMin, double yMax, int zBins, double zMin, double zMax){
  if(xBins == 0){
	std::cout << "X Bins value must be a non-zero integer" << std::endl;
	return;
  } else if(yBins == 0){
	std::cout << "Y Bins value must be a non-zero integer" << std::endl;
	return;
  } else if (zBins == 0){
	std::cout << "Z Bins value must be a non-zero integer" << std::endl;
	return;
  }
  std::string histogramRef = "h_" + hID;
  histoMap3D[histogramRef] = new TH3F(hID.c_str(),hTitle.c_str(),xBins,xMin,xMax,yBins,yMin,yMax,zBins,zMin,zMax);
}

void NpolHistoManager::FillHistograms(std::string hID, double value){
  std::string histogramRef = "h_"+hID;
  (histoMap[histogramRef])->Fill(value);
}

void NpolHistoManager::FillHistograms(std::string hID, double xValue, double yValue){
  std::string histogramRef = "h_"+hID;
  (histoMap2D[histogramRef])->Fill(xValue,yValue);
}

void NpolHistoManager::FillHistograms(std::string hID, double xValue, double yValue, double zValue){
  std::string histogramRef = "h_"+hID;
  (histoMap3D[histogramRef])->Fill(xValue,yValue,zValue);
}

void NpolHistoManager::WriteHistograms(){
  std::cout << "Writing Histograms to File!" << std::endl;
  std::map<std::string,TH1F *>::iterator histoIt2;
  for(histoIt2 = histoMap.begin(); histoIt2 != histoMap.end(); histoIt2++) {
	if(histoIt2->second == NULL) continue;
	histoIt2->second->Write();
  }

  std::map<std::string,TH2F *>::iterator histoIt3;
  for(histoIt3 = histoMap2D.begin(); histoIt3 != histoMap2D.end(); histoIt3++) {
	if(histoIt3->second == NULL) continue;
	histoIt3->second->Write();
  }

  std::map<std::string,TH3F *>::iterator histoIt4;
  for(histoIt4 = histoMap3D.begin(); histoIt4 != histoMap3D.end(); histoIt4++) {
	if(histoIt4->second == NULL) continue;
	histoIt4->second->Write();
  }

}

void NpolHistoManager::ClearHistograms(){
  histoMap.clear();
  histoMap2D.clear();
  histoMap3D.clear();
}

void NpolHistoManager::DeleteHistograms(){
  std::map<std::string,TH1F *>::iterator histoIt;
  for(histoIt = histoMap.begin(); histoIt != histoMap.end(); histoIt++) {
	delete histoIt->second;
  }
  std::map<std::string,TH2F *>::iterator histoIt2;
  for(histoIt2 = histoMap2D.begin(); histoIt2 != histoMap2D.end(); histoIt2++) {
	delete histoIt2->second;
  }
  std::map<std::string,TH3F *>::iterator histoIt3;
  for(histoIt3 = histoMap3D.begin(); histoIt3 != histoMap3D.end(); histoIt3++) {
	delete histoIt3->second;
  }
  
}

void NpolHistoManager::OpenFile(TString OutputFile){
  delete outFile2;
  outFile2 = new TFile(OutputFile,"RECREATE");
}

void NpolHistoManager::CloseFile() {
	outFile2->Close();
	delete outFile2;
	outFile2 = NULL;
}
