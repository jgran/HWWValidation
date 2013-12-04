#ifndef EVENTMAKER_H
#define EVENTMAKER_H

// user include files
//
// class decleration
//

class EventMaker : public edm::EDProducer {
public:
     explicit EventMaker (const edm::ParameterSet&);
     ~EventMaker();

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

     edm::InputTag dcsTag_;
	   std::string aliasprefix_;
     bool isData_;
     //std::string datasetName_;
     //std::string CMS2tag_;
};


#endif
