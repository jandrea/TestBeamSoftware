#ifndef __UTILITY__HH
#define __UTILITY__HH

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "stdint.h"
#include "DataFormats.h"
class TFile;
using std::string;

namespace Utility {
  
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiter);
  void getChannelMaskedHits( std::vector<int>& vec, const std::vector<int>& mch );
  void getChannelMaskedClusters( std::vector<tbeam::cluster*>& vec, const std::vector<int>& mch );
  void getChannelMaskedStubs( std::vector<tbeam::stub*>& vec, const std::vector<int>& mch );
 
  int readStubWord( std::map<std::string,std::vector<unsigned int> >& stubids, const uint32_t sWord );
  TH1* getHist1D(const char* hname);
  TH1* getHist1D(const string& hname);

  template <class T>
  bool fillHist1D(const char* hname, T value, double w=1.0) 
  {
    TH1* h = getHist1D(hname);
    if (!h) return false;
    h->Fill(value, w);
    return true;
  }
  template <class T>
  bool fillHist1D(const string& hname, T value, double w=1.0) {
    return fillHist1D(hname.c_str(), value, w);
  }
   template <class T>
  void fillHistofromVec( const std::vector<T>& vec, const char* h) {
    for( unsigned int i = 0; i<vec.size(); i++ ) {
      fillHist1D(h,vec.at(i));
    }
  }
  void fill2DHistofromVec( const std::vector<int>& vecC0, const std::vector<int>& vecC1,const char* h);

  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  TH2* getHist2D(const char* hname);
  TH2* getHist2D(const std::string& hname);
  template <class T1, class T2>
  bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, double w=1.0) {
    TH2* h = getHist2D(hname);
    if (!h) return false;
    h->Fill(xvalue, yvalue, w);
    return true;
  }
  template <class T1, class T2>
  bool fillHist2D(const std::string& hname, T1 xvalue, T2 yvalue, double w=1.0) {
    return fillHist2D(hname.c_str(), xvalue, yvalue, w);
  }

  // ---------------------------------------------
  // Convenience routine for filling Profile histograms
  // ---------------------------------------------
  TProfile* getHistProfile(const char* hname);
  TProfile* getHistProfile(const std::string& hname);
  template <class T1, class T2>
  bool fillHistProfile(const char* hname, T1 xvalue, T2 yvalue) {
    TProfile* h = getHistProfile(hname);
    if (!h) return false;
    h->Fill(xvalue, yvalue);
    return true;
  }
  template <class T1, class T2>
  bool fillHistProfile(const std::string& hname, T1 xvalue, T2 yvalue) {
    return fillHistProfile(hname.c_str(), xvalue, yvalue);
  }
  
  //  ------------------------------------------ //
  //  Convenience routines for fitting residuals //
  //  ------------------------------------------ //

  static double FuncPol1Gaus(Double_t * x, Double_t * par){
    double xx = x[0];
    double norm = par[0];
    double mean = par[1];
    double sigma = par[2];
    double cte = par[3];
    double slope = par[4];
    double f = norm*exp(-0.5*((xx-mean)/sigma)*((xx-mean)/sigma));
    f = f + cte + slope * xx;
    return f;
  }

  static double funcStepGaus(Double_t * x, Double_t * par){
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

  static double FuncStepGausShift(Double_t * x, Double_t * par){
    double xx = x[0];
    double pitch = par[0];
    double sigma = par[1];
    double norm = par[2];
    double cte = par[3];
    double shift = par[4];
    double f =0;
    if (xx<shift) f = norm*(1+TMath::Erf((xx+pitch/2.-shift)/(sqrt(2)*sigma)));
    if (xx>shift) f = norm*(1-TMath::Erf((xx-pitch/2.-shift)/(sqrt(2)*sigma)));
    f += cte;
    return f;
  }
  static double normBackgroundFunction(Double_t* x , Double_t* par)
  {
    double xx = x[0];
    double b0 = par[0]; 
    double x0 = par[2];
    double sigma = par[3]; 
    double norm = par[1]/(std::sqrt(2*TMath::Pi())*sigma); 
    
    double f = b0 +  norm*TMath::Exp( -1*std::pow(xx-x0,2.0)/(2*std::pow(sigma,2)) );
    return f; 
  }
  
  static double signalFunction(Double_t* x , Double_t* par)
  {
    double Ntracks = par[0];
    double pitch = par[1];
    double offset = par[2];
    double sigma0 = par[3];
    double xx = x[0] +  offset;
    double f0 = (Ntracks/pitch)*0.5*( TMath::Erf((xx + 0.5*pitch)/sigma0) -  TMath::Erf((xx-0.5*pitch)/sigma0) );
    return f0;    
  }
  static double residualFeI4Function(Double_t* x , Double_t* par)
  {
    double b0 = par[0];
    double b1 = par[1];
    double a0 = par[2];
    double xOffset_b = par[3];
    double xOffset_s = par[4];
    double sigma_b = par[5];
    double sigma_s = par[6];
    double pitch = par[7];
    double xx_s = x[0] +  xOffset_s;
    double xx = x[0];
    double f0 = (b0/pitch)*0.5*( TMath::Erf((xx_s + 0.5*pitch)/sigma_s) -  TMath::Erf((xx_s-0.5*pitch)/sigma_s) );
    double f1 = a0 +  (b1/(std::sqrt(2*TMath::Pi())*sigma_b))*TMath::Exp( -1*std::pow(xx-xOffset_b,2.0)/(2*std::pow(sigma_b,2)) );
    double f = f0 + f1; 
    return f;
  }
  static double residualFeI4_total(Double_t* x , Double_t* par)
  {
    double N_norm = par[0];
    double b0 = par[1];
    double b1 = par[2];
    double a0 = (1-(b0+b1));
    double xOffset_b = par[3];
    double xOffset_s = par[4];
    double sigma_b = par[5];
    double sigma_s0 = par[6];
    double sigma_s1 = par[7];
    double pitch = par[8];
    double xx_s = x[0] +  xOffset_s;
    double xx = x[0];
    double pars_Background[4] = { 0.0 , b0*N_norm , xOffset_b , sigma_b};
    double pars_s0[4] = { N_norm*a0 , pitch , xOffset_s , sigma_s0};
    double pars_s1[4] = { N_norm*b1 , pitch , xOffset_s , sigma_s1};
    double f_s0 = signalFunction( x, pars_s0);
    double f_s1 = signalFunction( x, pars_s1);
    double f_b = normBackgroundFunction( x, pars_Background);
    return f_s0 +  f_s1 + f_b;
  }
  static double stepFunction(Double_t* x , Double_t* par)
  {
    double norm = par[0]; 
    double widthStep = par[1];
    double smear = 1.0/par[2];
    double offset = par[3];
    double xx = x[0] +  offset;
    double f = norm*( TMath::Erf((xx + widthStep)*smear) -  TMath::Erf((xx-widthStep)*smear) );
    return f;
  }
  static double backgroundFunction(Double_t* x , Double_t* par)
  {
    double xx = x[0];
    double f = par[0] + par[1]*TMath::Gaus(xx, par[2] , par[3]);
    return f; 
  }
  static double residualFunction(Double_t* x , Double_t* par)
  {
    double f_Signal = stepFunction(x , &par[0]);
    double f_Background = backgroundFunction( x , &par[4]);
    double f = f_Signal + f_Background ;
    return f;
  }
  static double correlationFunction(Double_t* x , Double_t* par)
  {
    return x[0]*par[1] + par[0];
  }
  //--------------------------------------------------------------------//

  void removeTrackDuplicates(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<double> *xTkNoOverlap, std::vector<double> *yTkNoOverlap);
  void removeTrackDuplicates(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<double> *slopeTk, std::vector<double> *xTkNoOverlap, std::vector<double> *yTkNoOverlap, std::vector<double> *slopeTkNoOverlap);
  void removeTrackDuplicates(const tbeam::TelescopeEvent *telEv, std::vector<tbeam::Track>& tkNoOverlap);
  
  void cutTrackFei4Residuals(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<int> *colFei4, std::vector<int> *rowFei4, std::vector<double> *xSelectedTk, std::vector<double> *ySelectedTk, double xResMean, double yResMean, double xResPitch, double yResPitch);
  void cutTrackFei4Residuals(std::vector<double> *xTk, std::vector<double> *yTk, std::vector<double> *slopeTk, std::vector<int> *colFei4, std::vector<int> *rowFei4, std::vector<double> *xSelectedTk, std::vector<double> *ySelectedTk, std::vector<double> *slopeSelectedTk, double xResMean, double yResMean, double xResPitch, double yResPitch);
  
  void cutTrackFei4Residuals(const tbeam::FeIFourEvent* fei4ev ,const std::vector<tbeam::Track>& tkNoOverlap, std::vector<tbeam::Track>& selectedTk, double xResMean, double yResMean, double xResPitch, double yResPitch, bool doClosestTrack);

  double extrapolateTrackAtDUTwithAngles(const tbeam::Track& track, double FEI4_z, double offset, double zPlane, double theta);
  std::pair<double, double> extrapolateTrackAtDUTwithAngles(const tbeam::Track& track, double FEI4_z, double offset_d0, double zDUT_d0, double deltaZ, double theta);

}
#endif
