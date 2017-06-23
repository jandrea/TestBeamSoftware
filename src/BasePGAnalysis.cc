/*!
        \file                BasePGAnalysis.cc
        \brief               Example user code for Analysis
        \author              Suvankar Roy Chowdhury, Rajarshi Bhattacharya
        \date                05/07/16
        Support :            mail to : suvankar.roy.chowdhury@cern.ch, rajarshi.bhattacharya@cern.ch
*/
#include "TROOT.h"
#include "TInterpreter.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include <map>
#include <utility>
#include <vector>
#include <sstream>
#include <fstream>

#include "BasePGAnalysis.h"
using std::vector;
using std::map;
BasePGAnalysis::BasePGAnalysis() :
  BeamAnaBase::BeamAnaBase()
{
}
void BasePGAnalysis::bookHistograms() {
  BeamAnaBase::bookHistograms();
  hist_->bookTrackMatchHistograms();
}

unsigned long int BasePGAnalysis::readEventsLimit(Int_t target) {

  Int_t run;
  unsigned long int max;
  ifstream _f;
  Bool_t check=false;
  //  _f.open("/afs/cern.ch/user/r/rossia/workspace/CMSSW_7_3_0_pre1/BeamBasePGAna/TDR_Code/TestBeamSoftware/interface/eventsLimit.dat");
  _f.open("interface/eventsLimit.dat");
  while(_f.good() && !check){
    _f>>run>>max;
    if(run==target) check=true;
    else max=10000000;
  }
  _f.close();
  return max;
}


void BasePGAnalysis::beginJob() {
  BeamAnaBase::beginJob();
  nEntries_ = analysisTree()->GetEntries();
  hist_ = outFile();
  setAddresses();
  bookHistograms();
  analysisTree()->GetEntry(0);
  getCbcConfig(condEv()->cwd, condEv()->window);
}
 
void BasePGAnalysis::eventLoop()
{
  std::string filename=inFile();
  run_ = std::stoi(filename.substr(filename.find(".root")-3,3));
  cout<<"Run Number : "<<run_<<endl;
  maxEvent = getMaxEvt();
  if(nEntries_<maxEvent || maxEvent==0) maxEvent=nEntries_;

   Long64_t nbytes = 0, nb = 0;
   cout << "#Events=" << nEntries_ <<" -->  MAX EVENTS TO BE PROCESSED : "<<maxEvent<<endl;
   hist_->fillHist1D("EventInfo","nevents", nEntries_);

   std::cout << "CBC configuration:: SW=" << stubWindow()
             << "\tCWD=" << cbcClusterWidth()
             << "\tOffset1="<< cbcOffset1() 
             << "\tOffset2" << cbcOffset2()
   << std::endl;
   Float_t elogangle;
   long int trkFid = 0;
   long int det0clsMatch = 0;
   long int det1clsMatch = 0;
   long int clsMatchboth = 0;
   long int clsMatchany = 0;
   long int recostubMatchD1 = 0;
   long int cbcstubMatchD1 = 0;
   long int cbcstubWrongD1 = 0;
   unsigned long int lastBadevent = 0; 
   int nMatchedCluster = 0;

   long int runSub = -1;
   int trkChip;


   for (Long64_t jentry=0; jentry<nEntries_ && jentry<maxEvent;jentry++) {
     clearEvent();
     Long64_t ientry = analysisTree()->GetEntry(jentry);
     if (ientry < 0) break;
     if (jentry%1000 == 0) {
       cout << " Events processed. " << std::setw(8) << jentry 
	    << endl;
       runSub++;
     }
     if(jentry==0) {
       if(condEv()->DUTangle>=999000) elogangle=(999000.-(condEv()->DUTangle))/10;
       hist_->fillHist1D("EventInfo","hvSettings", condEv()->HVsettings);
       hist_->fillHist1D("EventInfo","dutAngle", elogangle);
       hist_->fillHist1D("EventInfo","alignAngle", dutangle()*180/TMath::Pi());
       hist_->fillHist1D("EventInfo","vcth", condEv()->vcth);
       hist_->fillHist1D("EventInfo","offset", cbcOffset1());
       hist_->fillHist1D("EventInfo","offset", cbcOffset2());
       hist_->fillHist1D("EventInfo","window", stubWindow());
       hist_->fillHist1D("EventInfo","tilt", static_cast<unsigned long int>(condEv()->tilt));
       cout << "Alignment Parameters" << aLparameteres();
     }
     hist_->fillHist1D("EventInfo","isPeriodic",isPeriodic());
     hist_->fillHist1D("EventInfo","isGoodFlag",isGoodEvent());

     if(!isGoodEvent())   {
      lastBadevent = jentry; 
      continue;
     }

     if(doTelMatching() && fei4Ev()->nPixHits != 1)   continue;
     
     hist_->fillHist1D("EventInfo","condData", condEv()->condData);
     hist_->fillHist1D("EventInfo","tdcPhase", static_cast<unsigned int>(condEv()->tdcPhase));
      
      setDetChannelVectors();
      const auto& d0c0 = *det0C0();
      const auto& d0c1 = *det0C1();
      const auto& d1c0 = *det1C0();
      const auto& d1c1 = *det1C1();
      fillCommonHistograms();
      //Telescope Matching
      if(doTelMatching() && hasTelescope()) {
        hist_->fillHist1D("TrackMatch","nTrackParams",telEv()->nTrackParams);
        
        hist_->fillHist1D("TrackMatch", "trkcluseff", 0);
        //Residual Calculation Now moved to AlignmentAnalysis
        //std::vector<double>  xtkDet0, xtkDet1;
        //getExtrapolatedTracks(xtkDet0, xtkDet1);
        std::vector<tbeam::Track>  fidTrkcoll;
        getExtrapolatedTracks(fidTrkcoll);
        //hist_->fillHist1D("TrackMatch", "nTrackParamsNodupl", xtkDet0.size());
        hist_->fillHist1D("TrackMatch", "nTrackParamsNodupl", fidTrkcoll.size());
        if(fidTrkcoll.empty())    continue;
        bool trkClsmatchD0 = false;
        bool trkClsmatchD1 = false;
        bool smatchD1 = false;
        bool cbcmatchD1 = false;
        double minclsresD0 = 9999.;
        double minclsresD1 = 9999.;
        double minclsposD0 = 9999.;
        double minclsposD1 = 9999.;
        double minStubresC0 = 9999.;
        double minStubposC0 = 9999.;
        int minStubStripC0 = 999.;
        int minClusStripD0 = 999;
        int minClusStripD1 = 999;
        int minClusWD0 = 999;
        int minClusWD1 = 999;
 
        int minHitStripD0 = 999;
        int minHitStripD1 = 999;
        double minHitresStripD0 = 999.;
        double minHitresStripD1 = 999.;

        for(auto &tk : fidTrkcoll) {
          double x0 = tk.xtkDut0; 
          hist_->fillHist1D("TrackMatch","hposxTkDUT0",x0); 
          //matching at det0
          for(auto& h : d0c0) {
            double res = x0 - (h-nstrips()/2)*dutpitch();
            if(std::fabs(res) < std::fabs(minHitresStripD0)) {
              minHitresStripD0 = res; 
              minHitStripD0 = h;
            }
          }   
          for(auto& cl : dutRecoClmap()->at("det0C0") ) {
            double res = x0 - (cl.x-nstrips()/2)*dutpitch();
            if(std::fabs(res) <= 4*resDUT())   trkClsmatchD0 = true; 
            if(std::fabs(res) < std::fabs(minclsresD0))  {
              minclsresD0 = res;
              minclsposD0 = (cl.x-nstrips()/2)*dutpitch();
              minClusStripD0 = cl.x;
              minClusWD0 = cl.size;
            }
          }
          hist_->fillHist1D("TrackMatch","hminposClsDUT0",minclsposD0);
          hist_->fillHist1D("TrackMatch","minresidualDUT0_1trkfid", minclsresD0);
          hist_->fillHist1D("TrackMatch","clswidthDUT0_1trkfid", minClusWD0);
          hist_->fillHist2D("TrackMatch","minclsTrkPoscorrD0", x0/dutpitch() + nstrips()/2 , minClusStripD0);
          //matching at det1
          double x1 = tk.xtkDut1;
          hist_->fillHist1D("TrackMatch","hposxTkDUT1",x1); 
          for(auto& h : d1c0) {
            double res = x1 - (h-nstrips()/2)*dutpitch();
            if(std::fabs(res) < std::fabs(minHitresStripD1)) {
              minHitresStripD1 = res; 
              minHitStripD1 = h;
            }
          }   
          for(auto& cl : dutRecoClmap()->at("det1C0") ) {
            double res = x1 - (cl.x-nstrips()/2)*dutpitch();
            if(std::fabs(x1 - (cl.x-nstrips()/2)*dutpitch()) <= 4*resDUT())   trkClsmatchD1 = true;
            if(std::fabs(res) < std::fabs(minclsresD1))  {
              minclsresD1 = res;
              minclsposD1 = (cl.x-nstrips()/2)*dutpitch();
              minClusStripD1 = cl.x;
              minClusWD1 = cl.size;
            }
          }

          for(auto& s : dutRecoStubmap()->at("C0"))  {
            double sposres = x1 - (s.x-nstrips()/2)*dutpitch();
            if(std::fabs(sposres) <= 4*resDUT())  smatchD1 = true;  
            if(std::fabs(sposres) < std::fabs(minStubresC0)) {
              minStubresC0 = sposres;
              minStubposC0 = (s.x-nstrips()/2)*dutpitch();
              minStubStripC0 = s.x;
            }
          }
	  for(auto& scbc: cbcstubChipids()->at("C0")) {
	    trkChip =(int) (x1/dutpitch()+nstrips()/2)/128;
	    //	    cout<<(x1/dutpitch()+nstrips()/2)<<" "<<trkChip<<endl;
            if(trkChip == scbc) cbcmatchD1=true;
	    //else cout << (x1/dutpitch()+nstrips()/2) << " - CBC Chip " << scbc <<endl;
	  }


          hist_->fillHist1D("TrackMatch","hminposClsDUT1",minclsposD1);
          hist_->fillHist1D("TrackMatch","minresidualDUT1_1trkfid", minclsresD1);
          hist_->fillHist1D("TrackMatch","clswidthDUT1_1trkfid", minClusWD1);
          hist_->fillHist2D("TrackMatch","minclsTrkPoscorrD1", x1/dutpitch() + nstrips()/2, minClusStripD1);
          //for stub
          hist_->fillHist1D("TrackMatch","sminresidualC0_1trkfid", minStubresC0);
          hist_->fillHist1D("TrackMatch","hminposStub",minStubposC0);
          hist_->fillHist2D("TrackMatch","minstubTrkPoscorrD1_all", x1/dutpitch() + nstrips()/2, minStubStripC0);
          if(smatchD1)  hist_->fillHist2D("TrackMatch","minstubTrkPoscorrD1_matched", x1/dutpitch() + nstrips()/2, minStubStripC0);  
	  if(nStubscbcSword()>=1) hist_->fillHist2D("TrackMatch","minstubTrkPoscorrD1_CBC", x1/dutpitch() + nstrips()/2, minStubStripC0);  
       }

	int trkchip=-1;
	if( (int)((fidTrkcoll.at(0).xtkDut1/dutpitch()+nstrips()/2)/128) == 0 ){
	  hist_->fillHist1D("TrackMatch", "trkCBCeff", 0);
	  hist_->fillHist1D("TrackMatch","effChip0Vtdc_den",static_cast<unsigned int>(condEv()->tdcPhase));
	  trkchip=0;
	}
        if((int)((fidTrkcoll.at(0).xtkDut1/dutpitch()+nstrips()/2)/128) == 1 ){
	  hist_->fillHist1D("TrackMatch", "trkCBCeff", 2);
	  hist_->fillHist1D("TrackMatch","effChip1Vtdc_den",static_cast<unsigned int>(condEv()->tdcPhase));
	  trkchip=1;
	}
	if(cbcstubChipids()->at("C0").size()==1){
	  if(cbcstubChipids()->at("C0").at(0) == 0 && trkchip==0){
	    hist_->fillHist1D("TrackMatch", "trkCBCeff", 1);
	    hist_->fillHist1D("TrackMatch","effChip0Vtdc_num",static_cast<unsigned int>(condEv()->tdcPhase));
	  }
	  if(cbcstubChipids()->at("C0").at(0) == 1  && trkchip==1){
	    hist_->fillHist1D("TrackMatch", "trkCBCeff", 3);
	    hist_->fillHist1D("TrackMatch","effChip1Vtdc_num",static_cast<unsigned int>(condEv()->tdcPhase));
	  }
	  if(cbcstubChipids()->at("C0").at(0) == 1  && trkchip==0){
	    hist_->fillHist1D("TrackMatch", "trkCBCeff", 4);
	  }
	  if(cbcstubChipids()->at("C0").at(0) == 0  && trkchip==1){
	    hist_->fillHist1D("TrackMatch", "trkCBCeff", 5);
	  }
	}
	if(cbcstubChipids()->at("C0").size()>1) hist_->fillHist1D("TrackMatch", "trkCBCeff", 6);

        hist_->fillHist1D("TrackMatch", "trkcluseff", 3);
	trkFid++;
	hist_->fillHist1D("TrackMatch", "trkFidStability", runSub);
        hist_->fillHist1D("TrackMatch","effVtdc_den",static_cast<unsigned int>(condEv()->tdcPhase));
        if(trkClsmatchD0)   {
          det0clsMatch++;
          hist_->fillHist1D("TrackMatch", "trkcluseff", 4);
          hist_->fillHist1D("TrackMatch", "clsD0Stability", runSub);
        }
        if(trkClsmatchD1)   {
          det1clsMatch++;
          hist_->fillHist1D("TrackMatch", "trkcluseff", 5);
          hist_->fillHist1D("TrackMatch", "clsD1Stability", runSub);
        }
        if(trkClsmatchD0 || trkClsmatchD1)   clsMatchany++;
        if(trkClsmatchD0 && trkClsmatchD1)   {
          clsMatchboth++;
          hist_->fillHist1D("TrackMatch", "trkcluseff", 6);
          hist_->fillHist1D("TrackMatch", "clsBothStability", runSub);
          hist_->fillHist1D("TrackMatch","effVtdc_num",static_cast<unsigned int>(condEv()->tdcPhase));
        }
        if(smatchD1) {
          recostubMatchD1++;
          hist_->fillHist1D("TrackMatch", "trkcluseff", 8);
          hist_->fillHist1D("TrackMatch", "stubStability", runSub);
          hist_->fillHist1D("TrackMatch","effRECOVtdc_num",static_cast<unsigned int>(condEv()->tdcPhase));
        }
        if(!trkClsmatchD0 && !trkClsmatchD1)  {
          hist_->fillHist1D("TrackMatch", "trkcluseff", 7);
        }
	if(nStubscbcSword()==1 && cbcmatchD1){
	  cbcstubMatchD1++;
	  hist_->fillHist1D("TrackMatch", "trkcluseff", 9);
          hist_->fillHist1D("TrackMatch", "cbcStability", runSub);
          hist_->fillHist1D("TrackMatch","effCBCVtdc_num",static_cast<unsigned int>(condEv()->tdcPhase));
	}
	else if(nStubscbcSword()==1 && !cbcmatchD1){
	  cbcstubWrongD1++;
	  hist_->fillHist1D("TrackMatch", "trkcluseff", 10);
          hist_->fillHist1D("TrackMatch", "cbcStability", runSub);
          hist_->fillHist1D("TrackMatch","effCBCVtdc_num",static_cast<unsigned int>(condEv()->tdcPhase));
	}
	else if(nStubscbcSword()>1)  hist_->fillHist1D("TrackMatch", "trkcluseff", 11);

	if(nStubscbcSword()==1 && dutRecoClmap()->at("det0C0").size()==1 && dutRecoClmap()->at("det1C0").size()==1){
	  hist_->fillHist1D("TrackMatch/CBCcheck", "EventN", jentry);
	  hist_->fillHist1D("TrackMatch/CBCcheck", "NClusterDUT0", dutRecoClmap()->at("det0C0").size());
	  hist_->fillHist1D("TrackMatch/CBCcheck", "NClusterDUT1", dutRecoClmap()->at("det1C0").size());
	  for(auto& cl : dutRecoClmap()->at("det0C0") ) {
	    hist_->fillHist1D("TrackMatch/CBCcheck", "ClusterPosDUT0", cl.x);
	  }
	  for(auto& cl : dutRecoClmap()->at("det1C0") ) {
	    hist_->fillHist1D("TrackMatch/CBCcheck", "ClusterPosDUT1", cl.x);
	  }
	  for(auto& s : dutRecoStubmap()->at("C0"))  {
	    hist_->fillHist1D("TrackMatch/CBCcheck", "RECOStubPosDUT1", s.x);
	    hist_->fillHist1D("TrackMatch/CBCcheck", "RECOStubResDUT1", fidTrkcoll.at(0).xtkDut1 - (s.x-nstrips()/2)*dutpitch());
	  }
	  for(auto& cl0 : dutRecoClmap()->at("det0C0") ) {
	    for(auto& cl1 : dutRecoClmap()->at("det1C0") ) {
	      hist_->fillHist1D("TrackMatch/CBCcheck", "ClusterDiff", abs(cl1.x-cl0.x));
	      hist_->fillHist2D("TrackMatch/CBCcheck", "ClusDifvsRes",cl1.x-cl0.x,fidTrkcoll.at(0).xtkDut1 - (cl1.x-nstrips()/2)*dutpitch());
	    }
	  }

	}
       	if(smatchD1){
	  hist_->fillHist1D("TrackMatch/RECOcheck", "EventN", jentry);
	  hist_->fillHist1D("TrackMatch/RECOcheck", "NClusterDUT0", dutRecoClmap()->at("det0C0").size());
	  hist_->fillHist1D("TrackMatch/RECOcheck", "NClusterDUT1", dutRecoClmap()->at("det1C0").size());
	  for(auto& cl : dutRecoClmap()->at("det0C0") ) {
	    hist_->fillHist1D("TrackMatch/RECOcheck", "ClusterPosDUT0", cl.x);
	  }
	  for(auto& cl : dutRecoClmap()->at("det1C0") ) {
	    hist_->fillHist1D("TrackMatch/RECOcheck", "ClusterPosDUT1", cl.x);
	  }
	  for(auto& s : dutRecoStubmap()->at("C0"))  {
	    hist_->fillHist1D("TrackMatch/RECOcheck", "RECOStubPosDUT1", s.x);
	    hist_->fillHist1D("TrackMatch/RECOcheck", "RECOStubResDUT1", fidTrkcoll.at(0).xtkDut1 - (s.x-nstrips()/2)*dutpitch());
	  }	  
	  for(auto& cl0 : dutRecoClmap()->at("det0C0") ) {
	    for(auto& cl1 : dutRecoClmap()->at("det1C0") ) {
	      hist_->fillHist1D("TrackMatch/RECOcheck", "ClusterDiff", abs(cl1.x-cl0.x));
	      hist_->fillHist2D("TrackMatch/RECOcheck", "ClusDifvsRes",cl1.x-cl0.x,fidTrkcoll.at(0).xtkDut1 - (cl1.x-nstrips()/2)*dutpitch());
	    }
	  }
	}
	if(smatchD1 && nStubscbcSword()==1){
	  hist_->fillHist1D("TrackMatch/Bothcheck", "EventN", jentry);
	  hist_->fillHist1D("TrackMatch/Bothcheck", "NClusterDUT0", dutRecoClmap()->at("det0C0").size());
	  hist_->fillHist1D("TrackMatch/Bothcheck", "NClusterDUT1", dutRecoClmap()->at("det1C0").size());
	  for(auto& cl : dutRecoClmap()->at("det0C0") ) {
	    hist_->fillHist1D("TrackMatch/Bothcheck", "ClusterPosDUT0", cl.x);
	  }
	  for(auto& cl : dutRecoClmap()->at("det1C0") ) {
	    hist_->fillHist1D("TrackMatch/Bothcheck", "ClusterPosDUT1", cl.x);
	  }
	  for(auto& s : dutRecoStubmap()->at("C0"))  {
	    hist_->fillHist1D("TrackMatch/Bothcheck", "RECOStubPosDUT1", s.x);
	    hist_->fillHist1D("TrackMatch/Bothcheck", "RECOStubResDUT1", fidTrkcoll.at(0).xtkDut1 - (s.x-nstrips()/2)*dutpitch());
	  }	  
	  for(auto& cl0 : dutRecoClmap()->at("det0C0") ) {
	    for(auto& cl1 : dutRecoClmap()->at("det1C0") ) {
	      hist_->fillHist1D("TrackMatch/Bothcheck", "ClusterDiff", abs(cl1.x-cl0.x));
	      hist_->fillHist2D("TrackMatch/Bothcheck", "ClusDifvsRes",cl1.x-cl0.x,fidTrkcoll.at(0).xtkDut1 - (cl1.x-nstrips()/2)*dutpitch());
	    }
	  }
	}
	if(!smatchD1 && nStubscbcSword()==1){
	  hist_->fillHist1D("TrackMatch/CBCOnlycheck", "EventN", jentry);
	  hist_->fillHist1D("TrackMatch/CBCOnlycheck", "NClusterDUT0", dutRecoClmap()->at("det0C0").size());
	  hist_->fillHist1D("TrackMatch/CBCOnlycheck", "NClusterDUT1", dutRecoClmap()->at("det1C0").size());
	  for(auto& cl : dutRecoClmap()->at("det0C0") ) {
	    hist_->fillHist1D("TrackMatch/CBCOnlycheck", "ClusterPosDUT0", cl.x);
	  }
	  for(auto& cl : dutRecoClmap()->at("det1C0") ) {
	    hist_->fillHist1D("TrackMatch/CBCOnlycheck", "ClusterPosDUT1", cl.x);
	  }
	  for(auto& s : dutRecoStubmap()->at("C0"))  {
	    hist_->fillHist1D("TrackMatch/CBCOnlycheck", "RECOStubPosDUT1", s.x);
	    hist_->fillHist1D("TrackMatch/CBCOnlycheck", "RECOStubResDUT1", fidTrkcoll.at(0).xtkDut1 - (s.x-nstrips()/2)*dutpitch());
	  }	  
	  for(auto& cl0 : dutRecoClmap()->at("det0C0") ) {
	    for(auto& cl1 : dutRecoClmap()->at("det1C0") ) {
	      hist_->fillHist1D("TrackMatch/CBCOnlycheck", "ClusterDiff", abs(cl1.x-cl0.x));
	      hist_->fillHist2D("TrackMatch/CBCOnlycheck", "ClusDifvsRes",cl1.x-cl0.x,fidTrkcoll.at(0).xtkDut1 - (cl1.x-nstrips()/2)*dutpitch());
	    }
	  }
	}
      }   
   }//event loop
   //error(1/N )sqrt( k(1 − k/N )).
   std::cout << "\n#events with 1 fid trk(both)=" << trkFid
             << "\n#events with atleast 1 matched cluster with 1 fid trk(D0)=" << det0clsMatch
             << "\n#events with atleast 1 matched cluster with 1 fid trk(D1)=" << det1clsMatch
             << "\n#events with atleast 1 matched cluster with 1 fid trk(any)=" << clsMatchany
             << "\n#events with atleast 1 matched cluster with 1 fid trk(both)=" << clsMatchboth
             << "\n#events with atleast 1 matched reco stub in D1=" << recostubMatchD1
             << "\n#events with atleast 1 matched CBC stub in D1=" << cbcstubMatchD1
             << "\n#events with atleast 1 wrong CBC stub in D1=" << cbcstubWrongD1
             << "\n#RECO Abs Stub Efficiency=" << double(recostubMatchD1)/double(trkFid) << "\tError=" << TMath::Sqrt(recostubMatchD1*(1.- double(recostubMatchD1)/double(trkFid) ))/double(trkFid)
             << "\n#CBC Abs Stub Efficiency=" << double(cbcstubMatchD1)/double(trkFid) << "\tError=" << TMath::Sqrt(cbcstubMatchD1*(1.- double(cbcstubMatchD1)/double(trkFid) ))/double(trkFid)
	     << std::endl;
}

void BasePGAnalysis::clearEvent() {
  BeamAnaBase::clearEvent();
}

void BasePGAnalysis::endJob() {
  BeamAnaBase::endJob();
  hist_->closeFile();
}

BasePGAnalysis::~BasePGAnalysis(){
  delete hist_;
}
