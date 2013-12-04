#ifndef MUONMAKER_H
#define MUONMAKER_H

//
// class declaration
//

class MuonMaker : public edm::EDProducer {
public:
     explicit MuonMaker (const edm::ParameterSet&);

private:
     virtual void beginJob() ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;
  double muonIsoValuePF(const reco::Muon& mu, const reco::Vertex& vtx, float coner, float minptn, float dzcut, int filterId);
  
      // ----------member data ---------------------------
  edm::InputTag muonsInputTag;
  edm::InputTag pfCandsInputTag;
  edm::InputTag vtxInputTag;
  edm::InputTag showerTag_;
  edm::InputTag beamSpot_tag_;
  edm::InputTag src_;
  std::string tevMuonsName;
  std::string aliasprefix_;
  std::string branchprefix_;
  edm::Handle<reco::PFCandidateCollection> pfCand_h;
  edm::Handle<reco::VertexCollection> vertexHandle;


};


#endif
