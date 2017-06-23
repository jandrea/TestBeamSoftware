#ifndef EventAlignment_h
#define EventAlignment_h

#include "BeamAnaBase.h"
//#include "CommonAnalysisTools.h"

#include "TString.h"
#include "Utility.h"
#include "DataFormats.h"
#include "Histogrammer.h"
#include <algorithm>
#include "TVirtualFitter.h"
#include "TF1.h"
#include "TFitResult.h"
#include <iostream>


using namespace std;

typedef std::pair<int, int> SingleGroup_Offset ;
typedef std::pair<int, int> FoundStub;
typedef std::pair<FoundStub, double > FoundStubInfo;


typedef std::pair<double, double> FEI4_Residual ;
typedef std::map<std::string , FEI4_Residual> FEI4_Residuals ; 


typedef std::pair<int,std::vector<int>> MatchedTracks;
typedef std::vector<MatchedTracks> AllMatchedTracks ; 
typedef std::pair<double,double> DutHits; 
typedef std::pair<double,DutHits> DutHitsTracking; 
typedef std::pair<int,std::vector<DutHitsTracking>> ExtrapolatedHits ;
typedef std::vector<ExtrapolatedHits> AllExtrapolatedHits ; 


class TH1;
class EventAlignment : public BeamAnaBase 
{
	public:
		EventAlignment(int pNumberEvents_perGroup , int pMaxOffset );
		~EventAlignment();
		void histoBooker();
		void beginJob();
		void eventLoop();
		void clearEvent();
		void endJob();
	 
	private :
		//std::string outFile_;
		Histogrammer* hist_;
		unsigned long int nEntries_;

		int fNumberEvents_perGroup;
		int fMaxEventOffset; 
		double fFeI4_p0 = 60.3; 
		double fTot_Max = 25;
		double fFeI4_p1 = 148 ;
		double fFeI4_p2 = 81341.0 ; 

		double fL1_sigma; 
		double fL1_mean;
	  	double fFeI4_Q_mpv;
	    double fFeI4_Q_sigma;


		double fFeI4_Correlation_y0;
		double fFeI4_Correlation_m;
		double fFeI4_Correlation_y0_unc;
		double fFeI4_Correlation_m_unc;


		double fTkPointing_x;
		double fTkPointing_y;

		double fFeI4_Window_x;
		double fFeI4_Window_y;
		double fBkground_FeI4_x;
		double fBkground_FeI4_y;


		FEI4_Residual fResidualX;
		FEI4_Residual fResidualY;

		double fDet0_ZOffset;
  		double fDet1_ZOffset;	
  		double fDet0_Offset;
		double fDet1_Offset;
		double fDet0_Resolution;
		double fDet1_Resolution;
		double fDet0_Resolution_Unc;
		double fDet1_Resolution_Unc;

		int fFidTracks_Det0=0;
		int fFidTracks_Det1=0;
		int fFidTracks_wHits_Det0;
		int fFidTracks_wHits_Det1;

		int fHits_Det0=0; 
		int fHits_Det1=0; 
		int fClusters_Det0=0;
		int fClusters_Det1=0;

		int fStubs =0 ;
		int fStubsMatchedToTracks=0;
		int fFidTracks;
		int fFidTracks_Single_Det0;
		int fFidTracks_Single_Det1;
		int fFidTracks_Single_wHits_Det0;
		int fFidTracks_Single_wHits_Det1;

		double fTrackContamination;
		double fTrackContamination_Err;
		double fTrackContamination_1;
		double fBackground;
		double fBackground_Err;
  


		double FeI4_Q_to_ToT(double Q){return fFeI4_p0*(fFeI4_p1+Q)/(fFeI4_p2+Q); };
		double FeI4_ToT_to_Q(double ToT){return (ToT*fFeI4_p2 - fFeI4_p0*fFeI4_p1)/(fFeI4_p0 - ToT); };

		void removeDuplicateTracks(std::vector<int>& idTks);

		std::vector<double> getRecoStubSeedWidths();
		std::vector<double> getRecoStubMatchWidths();
		std::vector<double> getRecoStubSeeds();
		std::vector<double> getRecoStubMatches();
		std::vector<double> getCbcStubPositions();
		std::vector<double> getRecoStubPositions();
		std::vector<double> getStubPositions();
		std::vector<double> getHitPositions( vector<int> detHits );
		std::vector<double> getClusterPositions( int pDetectorID  , int maxClusterWidth = 3);
		std::vector<int> getClusterWidths( int pDetectorID , int maxClusterWidth = 3);

		void lv1FeI4();
		tbeam::FeIFourEvent* sortFei4Hits(tbeam::FeIFourEvent* pFei4Ev );
		void clusterFeI4Hits( tbeam::FeIFourEvent* pFei4Ev , std::vector<double>& xPosClusters,std::vector<double>& yPosClusters, std::vector<double>& sizeClusters , std::vector<double>& totClusters , std::vector<double>& lv1Clusters , int minWidthCluster );
    	void CorrelateFeI4Tracks(int groupID, int evOffset );
		void CorrelateDetFeI4(int groupID, int evOffset );
		int MatchDetFeI4(int groupID, int prevOffset, double acceptanceWindow);
		int MatchDetFeI4_v2(int groupID, int prevOffset , double pCorrelation = 10 , double minCorrelation = 0.0 );
    
        double returnStrip(double xHit , double pitch = 0 );
        double returnHit( double iStrip, double pitch = 0 );
        void clusterCbcHits(std::vector<double> xHits, std::vector<double>& xPosClusters, std::vector<double>& sizeClusters, double maxClusterWidth=3);
        
        void sortStubs(std::vector<FoundStubInfo>& xStubs_Seeds , std::vector<FoundStubInfo>& xStubs_Matched , int correlationWindow_Strips = 5);
		uint8_t findStubs(std::vector<FoundStubInfo>& xStubs_Seeds , std::vector<FoundStubInfo>& xStubs_Matched , std::vector<double> xHits_det0 , std::vector<double> xHits_det1 , int correlationWindow_Strips=5 , bool halfStripRes = false);
        int findLongestCluster(std::vector<FoundStub> xStubs );
        
        void dut_RecoStubEfficiency( int pGroupID ,  int cWindow  =5 );
        void dutOnly_Efficiency(int pGroupID );
        void dutOnlyAnalysis(SingleGroup_Offset gOffset);
        void removeDuplicateHits( AllExtrapolatedHits& pExtrapolatedHits );

        // w/ tracking 
        void clearTrackMatchhistograms();
        void clearDUThistograms(const char* buffer);//modified from 
        void clearHistograms();
        bool checkHitIsFid( double xHit, int pDetectorID );
        double ExtrapolateFeI4Hit_DUT( int itrk ,int pDetectorID );
        bool checkHits(int pDetectorID , vector<int> detHits, double xTrk , int& pNumHits , std::vector<double>& residualsAbs , std::vector<double>& residuals );
        void CalculateHitResidual(AllExtrapolatedHits pExtrapolatedHits);
		void CalculateResolution(int pDetectorID);
        void CalculateEfficiency_Hits( AllExtrapolatedHits pExtrapolatedHits, double acceptanceWindow = 2.0*0.09 );
  		AllExtrapolatedHits MatchFeI4Hits(int groupID, int evOffset , double windowSize); 
		void alignX(std::vector<SingleGroup_Offset> gOffsets, double windowSize=1.0);
        AllExtrapolatedHits CorrelateFeI4Tracks(int groupID, int evOffset, double windowSizeX , double windowSizeY); 
        TF1* alignFeI4(TH1D* htmp , double pitch , double xExtra);
        TFitResultPtr fitFeI4Residual(TH1D* htmp ,  TH1D* htmp_bkg  , double pitch );
  		TFitResultPtr fitFeI4_Signal( TH1D* htmp , double pitch  );
  		TFitResultPtr fitFeI4_Background( TH1D* htmp , double pitch  );
  		void doFeI4TrackFit_HighIntensity();
  		
  		void estimateBackground(TFitResult* fitResult , double xMin , double xMax , double& meanIntegral , double& errorIntegral );
  		void estimateTrackContamination(double windowSize);
        void CalculateDUTefficiency( AllExtrapolatedHits pExtrapolatedHits, double acceptanceWindow = 0 , double correlationWindow  = 5);
  		void CalculateEfficiency(double matchingWindowSize=1.0);
        
		void FindDetFeI4CorrelationWindow(int nEvents=2000);
	public:
		std::vector<SingleGroup_Offset> alignEvents(std::string eAlignmentFile );
        private:
          tbeam::alignmentPars  alParsLocal_;
};
#endif
