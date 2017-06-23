#ifndef DataFormats_h
#define DataFormats_h
#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include<ostream>
#include<vector>
#include<map>
#include<string>
#include <string>
#include<stdint.h>
using std::vector;
using std::map;
using std::string;

namespace tbeam {
  class alignmentPars {
   public:
     alignmentPars();
     alignmentPars(const double fei4Z, const double resfei4X, const double resfei4Y, const double offxfei4, const double offyfei4, 
                   const double off_d0, const double zDUT_d0, const double dZ, const double t);
     ~alignmentPars(){}
     double residualSigmaFEI4x() const { return resSigmaFEI4X; } 
     void residualSigmaFEI4x(const double fx){ resSigmaFEI4X = fx; } 

     double residualSigmaFEI4y() const { return resSigmaFEI4Y; }
     void residualSigmaFEI4y(const double fy) { resSigmaFEI4Y = fy; }

     double offsetFEI4x() const { return offFEI4X; }
     void offsetFEI4x(const double ofx) { offFEI4X = ofx; }

     double offsetFEI4y() const { return offFEI4Y; }
     void offsetFEI4y(const double ofy) { offFEI4Y = ofy; }

     double FEI4z() const { return FEI4_z; }
     void FEI4z(const double fz) { FEI4_z = fz; }
 
     double d0Z() const { return d0_chi2_min_z; }
     void d0Z(const double dz) { d0_chi2_min_z = dz; }

     double d1Z() const { return d1_chi2_min_z; }
     void d1Z(const double dz) { d1_chi2_min_z = dz; }

     double d0Offset() const { return d0_Offset_aligned; } 
     void d0Offset(const double doff) { d0_Offset_aligned = doff; } 

     double d1Offset() const { return d1_Offset_aligned; } 
     void d1Offset(const double doff) { d1_Offset_aligned = doff; } 

     double deltaZ() const { return deltaz; }
     void deltaZ(const double dz) { deltaz = dz; }

     double theta() const { return angle;} 
     void theta(const double t) { angle = TMath::Pi()*t/180.; } 
    
     void setD1parametersfromD0();
     friend std::ostream& operator<< ( std::ostream& out, const tbeam::alignmentPars& ev );
   private:
     double FEI4_z;
     double d0_chi2_min_z;
     double d1_chi2_min_z;
     double d0_Offset_aligned;
     double d1_Offset_aligned;
     double deltaz;
     double angle;
     double resSigmaFEI4X;
     double resSigmaFEI4Y;
     double offFEI4X;
     double offFEI4Y;
  };
  class cluster;
  class stub;
  class cluster : public TObject {
   public:
     cluster();
     ~cluster();
     uint16_t x;
     float fx;
     uint16_t size;
     //std::vector<tbeam::stub *> stubs;
     ClassDef(cluster,1)
  };
  class stub : public TObject {
   public:
     stub();
     stub(const tbeam::stub& t);
     tbeam::cluster * seeding;  // Bottom sensor cluster
     tbeam::cluster * matched;  // Top sensor cluster
     uint16_t x;        // Position of the stub (bottom sensor)
     float fx;        // Position of the stub (bottom sensor) in 1/2 strip precision
     int16_t direction; // Direction of the stub (cl0-cl1)
     ClassDef(stub,1)
  };
  class cbc : public TObject {
   public:
     cbc();
     uint16_t pipelineAdd;
     uint8_t status;
     uint8_t error;
     ClassDef(cbc,1)
  };
  
  class dutEvent : public TObject {
   public:
     dutEvent();
     dutEvent(const tbeam::dutEvent& t);
     virtual ~dutEvent();
     std::map < std::string, std::vector <tbeam::cluster*> > clusters;
     std::map< std::string,std::vector<int> > dut_channel;
     std::map< std::string,std::vector<int> > dut_row;
     std::vector <tbeam::stub*> stubs;
     std::vector <tbeam::stub*> fstubs;//stubs with floating point precision of position
     uint32_t stubWord;
     uint32_t stubWordReco;
     //bool isGood;
     ClassDef(dutEvent,1)
  };
  class condEvent : public TObject {
   public:
    unsigned int run;
    unsigned int lumiSection;
    unsigned int event;
    unsigned long long time;
    unsigned long int unixtime;
    unsigned int tdcPhase;
    unsigned int HVsettings;
    unsigned int DUTangle;
    uint32_t window;
    uint32_t offset;
    uint32_t cwd;
    uint32_t tilt;
    uint32_t vcth;
    uint32_t stubLatency;
    uint32_t triggerLatency;
    int condData;
    int glibStatus;
    std::vector<tbeam::cbc> cbcs;
    condEvent();
    virtual ~condEvent(){}
    ClassDef(condEvent,1)
  };

  class TelescopeEvent : public TObject {
    public :
     Int_t           nTrackParams;
     Int_t           euEvt;
     std::vector<double>  *xPos;
     std::vector<double>  *yPos;
     std::vector<double>  *dxdz;
     std::vector<double>  *dydz;
     std::vector<int>     *trackNum;
     std::vector<int>     *iden;
     std::vector<double>  *chi2;
     std::vector<double>  *ndof;
     TelescopeEvent();
     TelescopeEvent(const TelescopeEvent& t);
     virtual ~TelescopeEvent();
     ClassDef(TelescopeEvent,1)
  };
  
  class FeIFourEvent : public TObject {
    public:
     Int_t                nPixHits;
     Int_t                euEvt;
     std::vector<int>     *col;
     std::vector<int>     *row;
     std::vector<int>     *tot;
     std::vector<int>     *lv1;
     std::vector<int>     *iden;
     std::vector<int>     *hitTime;
     std::vector<double>  *frameTime;
     FeIFourEvent();
     FeIFourEvent(const FeIFourEvent& f);
     virtual ~FeIFourEvent();
     ClassDef(FeIFourEvent,1)
  };
  
  class Track : public TObject {
    public :
      int trkIndex;
      double xPos;
      double yPos;
      double dxdz;
      double dydz;
      double chi2;
      double ndof;
      double xtkDut0;
      double xtkDut1;
      double ytkDut0;
      double ytkDut1;
      Track();
      Track(int i, double x, double y, double sx, double sy, double c2, double ndf);
      Track(int i, double x, double y, double sx, double sy, double c2, double ndf, double xtk0, double xtk1, double ytk0, double ytk1);
      Track(const tbeam::Track& t);
      virtual ~Track(){}
      ClassDef(Track,1)
  }; 

}
#endif
