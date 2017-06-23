#include <iostream>
#include <cstdlib>
#include <string>
#include "TROOT.h"
#include "TStopwatch.h"
#include "EventAlignment.h"
#include "argvparser.h"
using std::cout;
using std::cerr;
using std::endl;

using namespace CommandLineProcessing;

int main( int argc,char* argv[] )
{
   
  if(argc<2)  {
    std::cout << "Jobcard missing.\n./baselinReco <jobcardname>" << std::endl;
    return 1;
  }
  std::string jobfile = argv[1];
  //Let's roll
  TStopwatch timer;
  timer.Start();
  int cNevents_perGroup = 1000 ; 
  int cMaxOffset = 7;
  std::string eAlignmentFile = "";//FeI4_EventAlignment.dat";//"eAlignmentNov15.dat";//ealignmentNov15.txt";
  EventAlignment eAlign(cNevents_perGroup , cMaxOffset );
  std::cout << "Event Loop start" << std::endl;
  eAlign.readJob(jobfile);
  eAlign.beginJob();
  std::cout << "Event Loop start" << std::endl;
  // // return a vector of pairs <int,int> where the first element is the group id and the second is the event offset required 
  // // to optimize the correlation between the FeI4 clusters and DUT hits
  std::vector<SingleGroup_Offset> cGroupOffsets = eAlign.alignEvents(eAlignmentFile);
  eAlign.endJob();
  timer.Stop();
 
  cout << "Realtime/CpuTime = " << timer.RealTime() << "/" << timer.CpuTime() << endl;
  return 0;
}
