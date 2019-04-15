//***** Manager for Histograms and other Plots *****//

#ifndef Npol_Histogram_Manager_h
#define Npol_Hisogram_Manager_h

#include "TString.h"
#include <TROOT.h>
#include <TObject.h>
#include "TVector3.h"
#include "TVectorD.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <map>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

class NpolHistoManager {

public:
  static NpolHistoManager *GetInstance();
  void CreateHistograms(std::string hID, std::string hTitle,int xBins, double xMin, double xMax);
  void CreateHistograms(std::string hID, std::string hTitle,int xBins, double xMin, double xMax, int yBins, double yMin, double yMax);
  void CreateHistograms(std::string hID, std::string hTitle,int xBins, double xMin, double xMax, int yBins, double yMin, double yMax, int zBins, double zMin, double zMax);
  void ClearHistograms();
  void DeleteHistograms();
  void FillHistograms(std::string hID, double value);
  void FillHistograms(std::string hID, double xValue, double yValue);
  void FillHistograms(std::string hID, double xValue, double yValue, double zValue);
  void WriteHistograms();
  void OpenFile(TString OutputFile);
  void CloseFile();
  
private:
  NpolHistoManager();
  ~NpolHistoManager();
  
private:
  static NpolHistoManager *HistoMan;

public:
  TFile *outFile2 = NULL;
  std::map<std::string,TH1F *> histoMap;
  std::map<std::string,TH2F *> histoMap2D;
  std::map<std::string,TH3F *> histoMap3D;
  
};


#endif
