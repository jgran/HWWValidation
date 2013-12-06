#include "HWWValidation/HWWBase/interface/monitor.h"
#include "HWWValidation/HWWBase/interface/HWW.h"
#include "TH1F.h"
#include "TFile.h"
#include <fstream>
#include <string>

MonitorEventId::MonitorEventId(HWW& hww){
  run = hww.evt_run();
  event = hww.evt_event();
  lumi = hww.evt_lumiBlock();
}

MonitorEventId::MonitorEventId(){
  run = 0;
  event = 0;
  lumi = 0;
}

Entry::Entry()
{
  for (unsigned int i=0; i<5; ++i){
    nhyp[i] = 0;
    nevt[i] = 0;
    seen[i] = false;
  }
}

void hypo_monitor::count(HWW& hww, HypothesisType type, const char* name, double weight)
{
  std::vector<Entry>::iterator itr = counters.begin();
  while (itr != counters.end() && itr->name != name) itr++;
  Entry* entry(0);
  if ( itr == counters.end() ){
    counters.push_back(Entry());
    entry = &counters.back();
    entry->name = name;
  } else {
    entry = &*itr;
  }
  MonitorEventId id(hww);
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

void hypo_monitor::print() const
{
  std::cout << "Total number of processed events: \t" << nEvtProcessed << std::endl;
  ofstream out_file; out_file.open("cutflow.txt");
  for ( unsigned int i=0; i<counters.size(); ++i ){

    if( (counters[i].name).compare("njets == 0"     ) == 0) out_file << "-------------------------" << endl;
    if( (counters[i].name).compare("njets == 1"     ) == 0) out_file << "-------------------------" << endl;
    if( (counters[i].name).compare("njets == 2 or 3") == 0) out_file << "-------------------------" << endl;
    out_file << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		      counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) 	      
	      << std::endl;

    if( (counters[i].name).compare("njets == 0"     ) == 0) std::cout << "-------------------------" << endl;
    if( (counters[i].name).compare("njets == 1"     ) == 0) std::cout << "-------------------------" << endl;
    if( (counters[i].name).compare("njets == 2 or 3") == 0) std::cout << "-------------------------" << endl;
    std::cout << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		      counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) 	      
	      << std::endl;

    //store info about all the events passing each selection
    std::ofstream cut_file(Form("cut/cut-%d.txt",i));
    cut_file << Form("%-40s \tnevts: %u/%u/%u/%u/%u", counters[i].name.c_str(),
		     counters[i].nevt[MM],counters[i].nevt[EE],counters[i].nevt[EM],counters[i].nevt[ME],counters[i].nevt[ALL]) << "\n";

    for ( std::vector<MonitorEventId>::const_iterator id=counters[i].events.begin();
	  id!=counters[i].events.end(); ++id ){
      cut_file << id->run << "\t" << id->lumi << "\t" << id->event <<"\n";
    }
    cut_file.close();
  }
}

void hypo_monitor::makeHistograms() const
{
  TH1F* hist[4];
  TFile* outfile = new TFile("cutflow_hists.root", "RECREATE");
  outfile->cd();
  float denom = 0.0;
  float num = 0.0;
  for (unsigned int i=0; i<4; i++){
    hist[i]  = new TH1F(Form("cutflow_%s", HypothesisTypeName(i)), 
			Form("Relative Efficiency %s", HypothesisTypeName(i)), counters.size(), 0, counters.size() );	
    for (unsigned int j=0; j<counters.size(); ++j){
      if(j==0) denom = counters[0].nevt[i];//first cut will have efficiency of 1.0
      if(j>0)  denom = counters[j-1].nevt[i];//measure efficiency relative to previous cut
      num = counters[j].nevt[i]; 
      hist[i]->GetXaxis()->SetBinLabel(j+1,counters[j].name.c_str());
      if(denom==0) continue;
      hist[i]->SetBinContent(j+1,num/denom);
      float error = sqrt( (num/pow(denom,2))*(1 - num/denom) ); //binomial error
      hist[i]->SetBinError(j+1,error);
    }
    hist[i]->Write();
  }
}
