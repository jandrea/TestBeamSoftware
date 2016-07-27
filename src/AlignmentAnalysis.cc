/*!
 *         \file                AlignmentAnalysis.cc
 *         \brief               Alignment analysis, written starting from BaselineAnalysis
 *         \author              Nicolas Chanon
 *         \date                20/07/16
 *         Support :            mail to : nicolas.pierre.chanon@cern.ch
 * */

#include "TROOT.h"
#include "TInterpreter.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include <map>
#include <utility>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include "AlignmentAnalysis.h"

using std::vector;
using std::map;
using namespace std;

AlignmentAnalysis::AlignmentAnalysis(const std::string inFilename,const std::string outFilename) :
  BeamAnaBase::BeamAnaBase(),
  outFile_(outFilename),
  isProduction_(false)
{
  std::cout << "EDM Tuple Analyzer initialized with the following options\n" 
            << "Infile: " << inFilename
            << "\nOutFile: " << outFile_
            << std::endl; 
  if( setInputFile(inFilename) == 0 ) {
    std::cout << "Empty Chain!!";
    exit(1);
  }
  
  int pos1 = inFilename.find ("AnalysisTree_");
  pos1 += 13;
  int pos2 = inFilename.find (".root");
  runNumber = inFilename.substr(pos1, pos2-pos1);
  std::cout << "Run Number="<<runNumber<<endl;
  
  nEntries_ = analysisTree()->GetEntries();
  hist_ = new Histogrammer(outFile_);
  beginJob();
  
}
void AlignmentAnalysis::bookHistograms() {
  hist_->bookEventHistograms();
  hist_->bookTrackMatchingHisto();
}

void AlignmentAnalysis::beginJob() {
  setAddresses();
  bookHistograms();
  analysisTree()->GetEntry(0);
  getCbcConfig(condEv()->cwd, condEv()->window);
}

bool AlignmentAnalysis::setRunMode(const bool rMode) {
  isProduction_ = rMode;
}
void AlignmentAnalysis::eventLoop()
{
  if (!doTelMatching() || !hasTelescope()) return;
  
  float DUT_z = 435.0;
  float DUT_z_try = 400;
  float FEI4_z = 724.0;
  float pitch = 90.;
  int nMaxChannels = 254;
  
  Long64_t nbytes = 0, nb = 0;
  cout << "#Events=" << nEntries_ << endl;
  hist_->fillHist1D("EventInfo","nevents", nEntries_);
  
  std::cout << "CBC configuration:: SW=" << stubWindow()
	    << "\tCWD=" << cbcClusterWidth()
             << "\tOffset1="<< cbcOffset1() 
	    << "\tOffset2" << cbcOffset2()
	    << std::endl;
  //First Loop over events-inject z and compute residual
  //evaluate the best z alignment
  for (Long64_t jentry=0; jentry<nEntries_;jentry++) {
    clearEvent();
    Long64_t ientry = analysisTree()->GetEntry(jentry);
    if (ientry < 0) break;
    if (jentry%1000 == 0) {
       cout << " Events processed. " << std::setw(8) << jentry 
	    << endl;
    }
    if(jentry==0) {
      hist_->fillHist1D("EventInfo","hvSettings", condEv()->HVsettings);
      hist_->fillHist1D("EventInfo","dutAngle", condEv()->DUTangle);
      hist_->fillHist1D("EventInfo","vcth", condEv()->vcth);
      hist_->fillHist1D("EventInfo","offset", cbcOffset1());
      hist_->fillHist1D("EventInfo","offset", cbcOffset2());
      hist_->fillHist1D("EventInfo","window", stubWindow());
      hist_->fillHist1D("EventInfo","tilt", static_cast<unsigned long int>(condEv()->tilt));
    }
    hist_->fillHist1D("EventInfo","isPeriodic",isPeriodic());
    hist_->fillHist1D("EventInfo","isGoodFlag",isGoodEvent());
    
    if(!isGoodEvent())   continue;
    
    hist_->fillHist1D("EventInfo","condData", condEv()->condData);
    hist_->fillHist1D("EventInfo","tdcPhase", static_cast<unsigned int>(condEv()->tdcPhase));
    
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();      
    if (telEv()->nTrackParams==1){
      if (d0c0.size()==1){
	float xTkAtDUT = telEv()->xPos->at(0) + (DUT_z-FEI4_z)*telEv()->dxdz->at(0);
	float xDUT = (d0c0.at(0) - nMaxChannels/2) * pitch / 1000.;
	hist_->fillHist1D("TrackMatching","d0_1tk1Hit_diffX", xDUT-xTkAtDUT);
	DUT_z_try = 200; //300
	for (int iz=0; iz<50; iz++){
	  DUT_z_try += 10.;
	  //DUT_z_try += (float)(iz*5);
	  xTkAtDUT = telEv()->xPos->at(0) + (DUT_z_try-FEI4_z)*telEv()->dxdz->at(0);
	  hist_->fillHist1D("TrackMatching",Form("d0_1tk1Hit_diffX_iz%i", iz), xDUT-xTkAtDUT);
	}
      }
      if (d1c0.size()==1){
	float xTkAtDUT = telEv()->xPos->at(0) + (DUT_z-FEI4_z)*telEv()->dxdz->at(0);
	float xDUT = (d1c0.at(0) - nMaxChannels/2) * pitch / 1000.;
	hist_->fillHist1D("TrackMatching","d1_1tk1Hit_diffX", xDUT-xTkAtDUT);
	DUT_z_try = 200;//300;
	for (int iz=0; iz<50; iz++){
	  //DUT_z_try += (float)(iz*5);
	  DUT_z_try += 10.;
	  xTkAtDUT = telEv()->xPos->at(0) + (DUT_z_try-FEI4_z)*telEv()->dxdz->at(0);
	  hist_->fillHist1D("TrackMatching",Form("d1_1tk1Hit_diffX_iz%i", iz), xDUT-xTkAtDUT);
	}
      }
    }
  }//End of First Loop
  
  cout << "Computing chi2" << endl;
  float* sigma_d0 = new float[50];
  float* sigma_d1 = new float[50];
  float* offset_err_d0 = new float[50];
  float* offset_err_d1 = new float[50];
  
  std::pair<float, float> line_d0 = GetOffsetVsZ("d0", &sigma_d0, &offset_err_d0);
  std::pair<float, float> line_d1 = GetOffsetVsZ("d1", &sigma_d1, &offset_err_d1);
  
  cout << "d0, offset="<<line_d0.first <<" * z_DUT + " << line_d0.second << endl;
  cout << "d1, offset="<<line_d1.first <<" * z_DUT + " << line_d1.second << endl;
  
  float resTelescopeSquared = (90*90/12. + 3.5*3.5) / 1000. /1000.;
  float resTelescope = sqrt(0.090*0.090/12. + 0.0035*0.0035);
  float chi2_d0[50];
  float chi2_d1[50];
  float Nevent_d0_window[50];
  float Nevent_d1_window[50];
  for (int iz=0; iz<50; iz++) {
    chi2_d0[iz]=0;
    chi2_d1[iz]=0;
    Nevent_d0_window[iz]=0;
    Nevent_d1_window[iz]=0;
  }
  
  int nEv1tk1hit_d0=0;
  int nEv1tk1hit_d1=0;
   
  //2nd loop over events- computing chi2
  for (Long64_t jentry=0; jentry<nEntries_;jentry++) {
    clearEvent();
    Long64_t ientry = analysisTree()->GetEntry(jentry);
    if (ientry < 0) break;
    if (jentry%1000 == 0) {
      cout << " Events processed. " << std::setw(8) << jentry
	   << endl;
    }
    if(!isGoodEvent())   continue;
    if (telEv()->nTrackParams!=1) continue;
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();
    if (d0c0.size()==1){
      nEv1tk1hit_d0++;
      float xDUT = (d0c0.at(0) - nMaxChannels/2) * pitch / 1000.;
      for (int iz=0; iz<50; iz++){
	DUT_z_try = 200 + (float)(iz*10);
	float offset = line_d0.first * DUT_z_try + line_d0.second;
	float xTkAtDUT = telEv()->xPos->at(0) + (DUT_z_try-FEI4_z)*telEv()->dxdz->at(0) + offset;
	if (fabs(xTkAtDUT-xDUT) < 3*sigma_d0[iz]) {
	  chi2_d0[iz] += (xTkAtDUT-xDUT)/resTelescope * (xTkAtDUT-xDUT)/resTelescope;
	  Nevent_d0_window[iz] ++;
	}
      }
    }
    if (d1c0.size()==1){
      nEv1tk1hit_d1++;
      float xDUT = (d1c0.at(0) - nMaxChannels/2) * pitch / 1000.;
      for (int iz=0; iz<50; iz++){
	DUT_z_try = 200 + (float)(iz*10);
	float offset = line_d1.first * DUT_z_try + line_d1.second;
	float xTkAtDUT = telEv()->xPos->at(0) + (DUT_z_try-FEI4_z)*telEv()->dxdz->at(0) + offset;
	if (fabs(xTkAtDUT-xDUT) < 3*sigma_d1[iz]) {
	  chi2_d1[iz] += (xTkAtDUT-xDUT)/resTelescope * (xTkAtDUT-xDUT)/resTelescope;
	  Nevent_d1_window[iz]++;
	}
      }
    }
  }//end of 2nd loop over events
  
  
  for (int iz=0; iz<50; iz++){ 
    DUT_z_try = 200 + (float)(iz*10);
    chi2_d0[iz] /= Nevent_d0_window[iz];//nEv1tk1hit_d0;
    chi2_d1[iz] /= Nevent_d1_window[iz];//nEv1tk1hit_d1;
    cout << "z="<<DUT_z_try<<" chi2_d0="<<chi2_d0[iz]<<" chi2_d1="<<chi2_d1[iz]<<endl;
    if (chi2_d0[iz]>0) hist_->FillAlignmentOffsetVsZ("d0", "_chi2VsZ", iz, DUT_z_try, chi2_d0[iz], offset_err_d0[iz]/resTelescope);
    if (chi2_d1[iz]>0) hist_->FillAlignmentOffsetVsZ("d1", "_chi2VsZ", iz, DUT_z_try, chi2_d1[iz], offset_err_d1[iz]/resTelescope);
  }
  
  TF1* fParabolaChi2VsZ = new TF1("fParabolaChi2VsZ", "[0]*x*x+[1]*x+[2]", 200, 690);
  TH1* htmp = dynamic_cast<TH1F*>(hist_->GetHistoByName("TrackMatching","d0_chi2VsZ"));
  //TH1* htmp = (TH1F*) hist_->GetHistoByName("d0", "_chi2VsZ");
  htmp->Fit("fParabolaChi2VsZ");
  float d0_chi2_min = fParabolaChi2VsZ->GetMinimum();
  float d0_chi2_min_z = fParabolaChi2VsZ->GetMinimumX(); 
  cout << "d0 chi2 z="<<d0_chi2_min_z<<" chi2="<<d0_chi2_min<<endl;
  
  htmp = dynamic_cast<TH1F*>(hist_->GetHistoByName("TrackMatching","d1_chi2VsZ"));
  //htmp = (TH1F*) hist_->GetHistoByName("d1", "_chi2VsZ");
  htmp->Fit("fParabolaChi2VsZ");
  float d1_chi2_min = fParabolaChi2VsZ->GetMinimum();
  float d1_chi2_min_z = fParabolaChi2VsZ->GetMinimumX();
  cout << "d1 chi2 z="<<d1_chi2_min_z<<" chi2="<<d1_chi2_min<<endl;
  
  float DUT_z_aligned = DUT_z;
  float d0_Offset_aligned = line_d0.first * d0_chi2_min_z + line_d0.second;
  float d1_Offset_aligned = line_d1.first * d1_chi2_min_z + line_d1.second;
  
  ofstream fileAlignment;
  if(isProduction_) {
    fileAlignment.open("interface/alignmentMay16.txt", ios::out | ios::app);
  } else {
    std::string dumpFname = "alignmentMay16_" + runNumber + ".txt";
    fileAlignment.open(dumpFname.c_str(), ios::out);
  }
  if(fileAlignment) {
    fileAlignment << "Run="<< runNumber 
                  << ":zD0="<< d0_chi2_min_z 
                  << ":offsetD0=" << d0_Offset_aligned 
                  << ":zD1="<< d1_chi2_min_z 
                  << ":offsetD1=" << d1_Offset_aligned << endl;
    fileAlignment.close();
  } else std::cout << "Dump File could not be opened!!" << std::endl;

  //3rd loop over number of entries-get residual after alignment correction
  for (Long64_t jentry=0; jentry<nEntries_;jentry++) {
    clearEvent();
    Long64_t ientry = analysisTree()->GetEntry(jentry);
    if (ientry < 0) break;
    if (jentry%1000 == 0) {
      cout << " Events processed. " << std::setw(8) << jentry
	   << endl;
    }
    if(!isGoodEvent())   continue;
    if (telEv()->nTrackParams!=1) continue;
    setDetChannelVectors();
    const auto& d0c0 = *det0C0();
    const auto& d0c1 = *det0C1();
    const auto& d1c0 = *det1C0();
    const auto& d1c1 = *det1C1();
    
    if (d0c0.size()==1){
      float xDUT = (d0c0.at(0) - nMaxChannels/2) * pitch / 1000.;
      float xTkAtDUT = telEv()->xPos->at(0) + (d0_chi2_min_z-FEI4_z)*telEv()->dxdz->at(0) + d0_Offset_aligned;
      hist_->fillHist1D("TrackMatching","d0_1tk1Hit_diffX_aligned", xDUT-xTkAtDUT);
    }
    if (d1c0.size()==1){
      float xDUT = (d1c0.at(0) - nMaxChannels/2) * pitch / 1000.;
      float xTkAtDUT = telEv()->xPos->at(0) + (d1_chi2_min_z-FEI4_z)*telEv()->dxdz->at(0) + d1_Offset_aligned;
      hist_->fillHist1D("TrackMatching","d1_1tk1Hit_diffX_aligned", xDUT-xTkAtDUT);
    }
  }//end of 3rd loop over events

  //Fit Residuals Gaussian convulated with Step Function
  TF1* fGausResiduals = new TF1("fGausResiduals", "gaus", -10, 10);
  htmp = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatching","d0_1tk1Hit_diffX_aligned"));   
  //htmp = (TH1I*) hist_->GetHistoByName("d0", "_1tk1Hit_diffX_aligned");
  float center = ((float)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((float)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();//htmp->GetMean();
  htmp->SetAxisRange(center-0.3, center+0.3, "X");
  fGausResiduals->SetRange(center-0.3, center+0.3);
  fGausResiduals->SetParameter(1, center);
  fGausResiduals->SetParameter(2, 0.026);
  htmp->Fit("fGausResiduals");
  
  center = fGausResiduals->GetParameter(1);
  float rms = 5*fGausResiduals->GetParameter(2);
  htmp->SetAxisRange(center-rms, center+rms, "X");
  //fGausResiduals->SetRange(center-rms, center+rms);
  //htmp->Fit("fGausResiduals");
  
  TF1* fStepGaus = new TF1("FuncStepGaus", FuncStepGaus, -0.1, 0.1, 4);
  fStepGaus->SetLineWidth(2);
  fStepGaus->SetLineColor(kRed);
  fStepGaus->SetParLimits(0, 0., 0.1);
  fStepGaus->SetParameter(0, 0.09);
  fStepGaus->SetParLimits(1, 0.0, 0.05);
  fStepGaus->SetParameter(1, 0.01);
  fStepGaus->SetParameter(2, htmp->GetMaximum());
  fStepGaus->SetParLimits(3, 0, 10);
  fStepGaus->SetParameter(3, 0);
  
  htmp->Fit(fStepGaus);
  htmp = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatching","d1_1tk1Hit_diffX_aligned"));   
  //htmp = (TH1I*) hist_->GetHistoByName("d1", "_1tk1Hit_diffX_aligned");
  center = ((float)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((float)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();//htmp->GetMean();
  htmp->SetAxisRange(center-0.3, center+0.3, "X");
  fGausResiduals->SetRange(center-0.3, center+0.3);
  fGausResiduals->SetParameter(1, center);
  fGausResiduals->SetParameter(2, 0.026);
  htmp->Fit("fGausResiduals");
  
  center = fGausResiduals->GetParameter(1);
  rms = 5*fGausResiduals->GetParameter(2);
  htmp->SetAxisRange(center-rms, center+rms, "X");
  //fGausResiduals->SetRange(center-rms, center+rms);
  //htmp->Fit("fGausResiduals");
  fStepGaus->SetParameter(0, 0.09);
  fStepGaus->SetParameter(1, 0.01);
  fStepGaus->SetParameter(2, htmp->GetMaximum());
  fStepGaus->SetParameter(3, 0);
  htmp->Fit(fStepGaus);
}

std::pair<float, float> AlignmentAnalysis::GetOffsetVsZ(const char* det, float** sigma, float** offset_err){
  
  TF1* fGausExtractedX = new TF1("fGausExtractedX", "gaus", -10, 10);
  TH1* htmp;
  float offset[50];
  //float offset_err[50];
  float injectedZ[50];
  float center, rms;
  bool failedFit = false;
  
  for (int iz=0; iz<50; iz++){
    std::string hn = det + std::string("_1tk1Hit_diffX_iz") + std::to_string((iz));
    htmp = dynamic_cast<TH1I*>(hist_->GetHistoByName("TrackMatching", hn));   
    //htmp = (TH1I*) hist_->GetHistoByName(det, Form("_1tk1Hit_diffX_iz%i",iz));
    center = ((float)htmp->GetMaximumBin())*(htmp->GetXaxis()->GetXmax()-htmp->GetXaxis()->GetXmin())/((float)htmp->GetNbinsX()) + htmp->GetXaxis()->GetXmin();//htmp->GetMean();
    htmp->SetAxisRange(center-0.2, center+0.2, "X");
    fGausExtractedX->SetRange(center-0.2, center+0.2);
    fGausExtractedX->SetParameter(1, center);
    fGausExtractedX->SetParameter(2, 0.026);
    htmp->Fit("fGausExtractedX");
    
    center = fGausExtractedX->GetParameter(1);
    rms = 5*fGausExtractedX->GetParameter(2);
    (*sigma)[iz] = rms;
    htmp->SetAxisRange(center-rms, center+rms, "X");
    fGausExtractedX->SetRange(center-rms, center+rms);
    htmp->Fit("fGausExtractedX");
    
    offset[iz] = fGausExtractedX->GetParameter(1);
    (*offset_err)[iz] = fGausExtractedX->GetParError(1);
    injectedZ[iz] = 200 + (float)(iz*10);
    //failedFit = false;
    //if ((*offset_err)[iz]/offset[iz] > 1 || !((*offset_err)[iz]>0)) failedFit = true;
    //if (!failedFit) hist_->FillAlignmentOffsetVsZ(det, "_offsetVsZ", iz, injectedZ[iz], offset[iz], (*offset_err)[iz]);
  }
  
  for (int iz=0; iz<50; iz++){
    cout << "iz="<< iz<<" z="<<injectedZ[iz]<<" offset = "<<offset[iz]<< " +/- " <<(*offset_err)[iz] << endl;
    failedFit = false;
    if ((*offset_err)[iz]/offset[iz] > 1 || !((*offset_err)[iz]>0)) failedFit = true;
    if (!failedFit) hist_->FillAlignmentOffsetVsZ(det, "_offsetVsZ", iz, injectedZ[iz], offset[iz], (*offset_err)[iz]);
  }
  
  TF1* fOffsetVsZ = new TF1("fOffsetVsZ","[0]*x+[1]", injectedZ[0], injectedZ[49]);
  std::string hn = det + std::string("_offsetVsZ");
  htmp = dynamic_cast<TH1F*>(hist_->GetHistoByName("TrackMatching", hn));  
  //htmp = (TH1I*) hist_->GetHistoByName(det, "_offsetVsZ");
  htmp->Fit(fOffsetVsZ);
  htmp->SetAxisRange(htmp->GetMinimum(), htmp->GetMaximum(),"Y");
  
  std::pair<float, float> line;
  line.first = fOffsetVsZ->GetParameter(0);
  line.second = fOffsetVsZ->GetParameter(1);
  
  return line;
}

void AlignmentAnalysis::clearEvent() {
  BeamAnaBase::clearEvent();
}
void AlignmentAnalysis::endJob() {
  BeamAnaBase::endJob();
  hist_->closeFile();
}

AlignmentAnalysis::~AlignmentAnalysis(){
  delete hist_;
}