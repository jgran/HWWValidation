#ifndef PFMETMAKER_H
#define PFMETMAKER_H

//
// class decleration
//

class PFMETMaker : public edm::EDProducer {
public:
    explicit PFMETMaker (const edm::ParameterSet&);
    ~PFMETMaker();

private:
    virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------
    edm::InputTag pfMetInputTag;
    edm::InputTag pfMetCorInputTag;
	  std::string aliasprefix_;
    
};


#endif
