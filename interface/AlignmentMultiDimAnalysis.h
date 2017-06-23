#ifndef AlignmentMultiDimAnalysis_h
#define AlignmentMultiDimAnalysis_h

#include "BeamAnaBase.h"

#include "TString.h"
//#include "Utility.h"
#include "DataFormats.h"
#include "Histogrammer.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

class TH1;
class AlignmentMultiDimAnalysis : public BeamAnaBase {
 public:
  AlignmentMultiDimAnalysis();
  ~AlignmentMultiDimAnalysis();
  void beginJob();
  void eventLoop(); 
  void bookHistograms();
  void clearEvent();
  void endJob();
  double ComputeChi2(const double* x) const;
  double ComputeChi2BothPlanes(const double* x) const;

  static double FuncStepGaus(Double_t * x, Double_t * par){
    double xx = x[0];
    double pitch = par[0];
    double sigma = par[1];
    double norm = par[2];
    double cte = par[3];
    double f =0;
    if (xx<0) f = norm*(1+TMath::Erf((xx+pitch/2.)/(sqrt(2)*sigma)));
    if (xx>0) f = norm*(1-TMath::Erf((xx-pitch/2.)/(sqrt(2)*sigma)));
    f += cte;
    return f;
  }
  void doTelescopeAnalysis(tbeam::alignmentPars& aLp);
 private:
  std::string runNumber_;
  std::string outFile_;
  Histogrammer* hist_;
  int run_;
  unsigned long int nEntries_;
  unsigned long int maxEvent;
  bool isProduction_;
  std::string alignparFile_; 
  float zMin;
  float zStep;
  int zNsteps;

  bool doD0, doD1; 
  std::vector<tbeam::Track>  selectedTk_d0_1Hit;
  std::vector<tbeam::Track>  selectedTk_d1_1Hit;
  std::vector<tbeam::Track>  selectedTk_bothPlanes_1Cls; 
  std::vector<float> d0_DutXpos;
  std::vector<float> d1_DutXpos;
  std::vector<float> bothPlanes_DutXposD0;
  std::vector<float> bothPlanes_DutXposD1;
 
  ROOT::Math::Functor* toMinimize;
  ROOT::Math::Functor* toMinimizeBothPlanes;
  ROOT::Math::Functor* toMinimizeBothPlanesConstraint;
  ROOT::Minuit2::Minuit2Minimizer* minimizer;
  ROOT::Minuit2::Minuit2Minimizer* minimizerBothPlanes;
  ROOT::Minuit2::Minuit2Minimizer* minimizerBothPlanesConstraint;

  bool doConstrainDeltaOffset;
  tbeam::alignmentPars al;
};

#endif
