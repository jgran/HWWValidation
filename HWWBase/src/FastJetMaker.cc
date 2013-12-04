#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"

#include "HWWValidation/HWWBase/interface/FastJetMaker.h"


using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

FastJetMaker::FastJetMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
       branchprefix.replace(branchprefix.find("_"),1,"");

  produces<float>         (branchprefix+"rhoJEC"    ).setBranchAlias(aliasprefix_+"_rhoJEC"   );
  produces<float>         (branchprefix+"rho"       ).setBranchAlias(aliasprefix_+"_rho"      );

  // input tags
  rhoJEC_tag     = iConfig.getParameter<edm::InputTag>("rhoJEC_tag");
  rhoIso_tag     = iConfig.getParameter<edm::InputTag>("rhoIso_tag");
}


FastJetMaker::~FastJetMaker()
{
}

void  FastJetMaker::beginJob()
{
}

void FastJetMaker::endJob()
{
}

// ------------ method called to produce the data  ------------
void FastJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<float> evt_rho(new float);
  auto_ptr<float> evt_rhoIso(new float);
 
  edm::Handle<double> rhoH;
  iEvent.getByLabel( rhoJEC_tag , rhoH);

  edm::Handle<double> rhoIsoH;
  iEvent.getByLabel( rhoIso_tag , rhoIsoH);
  
  *evt_rho = *rhoH; 
  *evt_rhoIso = *rhoIsoH; 

  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos)
       branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(evt_rho            , branchprefix+"rhoJEC"     );
  iEvent.put(evt_rhoIso         , branchprefix+"rho"        );
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastJetMaker);
