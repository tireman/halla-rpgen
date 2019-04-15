//******* NPOL Pre-processing Class *******//


#ifndef Npol_Event_PreProcessing_h
#define Npol_Event_PreProcessing_h

#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"

class NpolEventPreProcessing {
  
public:
  NpolEventPreProcessing();
  ~NpolEventPreProcessing();
  static NpolEventPreProcessing *GetInstance();
  
  int GetAVNumber(const std::string &volName);
  int GetImprNumber(const std::string &volName);
  int GetPlacementNumber(const std::string &volName);
  void AnalyzerHitPosition(double hPos[],double lPos[], int detNums[]);
  void TaggerHitPosition(double hPos[],double lPos[], int detNums[]);
  void DeltaEarrayHitPosition(double hPos[],double lPos[], int detNums[]);
  void EarrayHitPosition(double hPos[],double lPos[], int detNums[]);
  void RotateNpolToG4(double hPos[], double RotAng);
  void RotateDetToNpol(double hPos[], int detNums[], double detAngle);
  void RotateG4ToRoot(double curPos[], double newPos[], double rotMat[][3]);
  static double NpolAng;
  
private:
  static NpolEventPreProcessing *PreProcess;
  TRandom3 *rand = new TRandom3();
};


#endif
