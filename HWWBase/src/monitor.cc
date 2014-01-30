#include <fstream>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "HWWValidation/HWWBase/interface/monitor.h"

//DQM
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


EventMonitor::MonitorEventId::MonitorEventId(){
  run = 0;
  event = 0;
  lumi = 0;
}

EventMonitor::Entry::Entry()
{
  for (unsigned int i=0; i<5; ++i){
    nhyp[i] = 0;
    nevt[i] = 0;
    nhyp_weighted[i] = 0.0;
    nevt_weighted[i] = 0.0;
    seen[i] = false;
  }
}

void EventMonitor::hypo_monitor::count(HypothesisType type, const char* name, double weight)
{
  std::vector<EventMonitor::Entry>::iterator itr = counters.begin();
  while (itr != counters.end() && itr->name != name) itr++;
  EventMonitor::Entry* entry(0);
  if ( itr == counters.end() ){
    counters.push_back(Entry());
    entry = &counters.back();
    entry->name = name;
  } else {
    entry = &*itr;
  }
  EventMonitor::MonitorEventId id;
  entry->nhyp[type]++;
  entry->nhyp[ALL]++;
  if (id != entry->lastEvent){
    if (keepEventList) entry->events.push_back(id);
    for (unsigned int i=0; i<5; ++i) entry->seen[i] = false;
    entry->nevt[type]++;
    entry->nevt[ALL]++;
    entry->nevt_weighted[type]+=weight;
    entry->nevt_weighted[ALL]+=weight;
    entry->seen[type] = true;
    entry->lastEvent = id;
  }
}

void EventMonitor::hypo_monitor::print() const
{
  std::cout << "Total number of processed events: \t" << nEvtProcessed << std::endl;
  ofstream out_file; out_file.open("cutflow.txt");
  for ( unsigned int i=0; i<counters.size(); ++i ){

    if( (counters[i].name).compare("njets == 0"     ) == 0) out_file << "-------------------------" << std::endl;
    if( (counters[i].name).compare("njets == 1"     ) == 0) out_file << "-------------------------" << std::endl;
    if( (counters[i].name).compare("njets == 2 or 3") == 0) out_file << "-------------------------" << std::endl;
    out_file << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		      counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) 	      
	      << std::endl;

    if( (counters[i].name).compare("njets == 0"     ) == 0) std::cout << "-------------------------" << std::endl;
    if( (counters[i].name).compare("njets == 1"     ) == 0) std::cout << "-------------------------" << std::endl;
    if( (counters[i].name).compare("njets == 2 or 3") == 0) std::cout << "-------------------------" << std::endl;
    std::cout << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		      counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) 	      
	      << std::endl;

    //store info about all the events passing each selection
    std::ofstream cut_file(Form("cut/cut-%d.txt",i));
    cut_file << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		     counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) << "\n";

    for ( std::vector<EventMonitor::MonitorEventId>::const_iterator id=counters[i].events.begin();
	  id!=counters[i].events.end(); ++id ){
      cut_file << id->run << "\t" << id->lumi << "\t" << id->event <<"\n";
    }
    cut_file.close();
  }
}
