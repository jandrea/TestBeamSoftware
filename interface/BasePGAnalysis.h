#ifndef BasePGAnalysis_h
#define BasePGAnalysis_h

#include "BeamAnaBase.h"

#include "TString.h"
//#include "Utility.h"
#include "DataFormats.h"
#include "Histogrammer.h"

class TH1;
class BasePGAnalysis : public BeamAnaBase {
 public:
  BasePGAnalysis();
  ~BasePGAnalysis();
  void beginJob();
  void eventLoop(); 
  void bookHistograms();
  void printEfficiency();
  void clearEvent();
  void endJob();
  unsigned long int readEventsLimit(Int_t target);
 private:
  //std::string outFile_;
  Histogrammer* hist_;
  int run_;
  unsigned long int nEntries_; 
  unsigned long int maxEvent; 
};
#endif
