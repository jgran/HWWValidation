#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "HWWValidation/HWWBase/interface/EventMaker.h"



typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventMaker::EventMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     produces<unsigned int>        (branchprefix+"run"            ).setBranchAlias(aliasprefix_+"_run"           );
     produces<unsigned int>        (branchprefix+"event"          ).setBranchAlias(aliasprefix_+"_event"         );
     produces<unsigned int>        (branchprefix+"lumiBlock"      ).setBranchAlias(aliasprefix_+"_lumiBlock"     );
     produces<int>                 (branchprefix+"bunchCrossing"  ).setBranchAlias(aliasprefix_+"_bunchCrossing" );
     produces<int>                 (branchprefix+"orbitNumber"    ).setBranchAlias(aliasprefix_+"_orbitNumber"   );
     produces<int>                 (branchprefix+"storeNumber"    ).setBranchAlias(aliasprefix_+"_storeNumber"   );
     produces<int>                 (branchprefix+"experimentType" ).setBranchAlias(aliasprefix_+"_experimentType");
     produces<float>               (branchprefix+"bField"         ).setBranchAlias(aliasprefix_+"_bField"        );
     produces<unsigned int>        (branchprefix+"detectorStatus" ).setBranchAlias(aliasprefix_+"_detectorStatus");
     produces<int>                 (branchprefix+"isRealData"     ).setBranchAlias(aliasprefix_+"_isRealData"    );
  
     dcsTag_ = iConfig.getParameter<edm::InputTag>("dcsTag");
     isData_ = iConfig.getParameter<bool>("isData") ;
}


EventMaker::~EventMaker() {}

void EventMaker::beginJob() {  
}

void EventMaker::endJob() {
}

// ------------ method called to produce the data  ------------
void EventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
     auto_ptr<unsigned int>         evt_run             (new unsigned int              );
     auto_ptr<unsigned int>         evt_event           (new unsigned int              );
     auto_ptr<unsigned int>         evt_lumiBlock       (new unsigned int              );
     auto_ptr<int>                  evt_bunchCrossing   (new int                       );
     auto_ptr<int>                  evt_orbitNumber     (new int                       );
     auto_ptr<int>                  evt_storeNumber     (new int                       );
     auto_ptr<int>                  evt_experimentType  (new int                       );
     auto_ptr<float>                evt_bField          (new float                     );
     auto_ptr<unsigned int>         evt_detectorStatus  (new unsigned int              );
     auto_ptr<int>                  evt_isRealData      (new int                       );
     
     *evt_run                       = iEvent.id().run()        ;
     *evt_event                     = iEvent.id().event()      ;
     *evt_lumiBlock                 = iEvent.luminosityBlock() ;
     *evt_bunchCrossing             = iEvent.bunchCrossing()   ;
     *evt_orbitNumber               = iEvent.orbitNumber()     ;
     *evt_storeNumber               = iEvent.eventAuxiliary().storeNumber();
     *evt_experimentType            = iEvent.experimentType()  ;
     *evt_isRealData                = iEvent.isRealData();

     edm::Handle<DcsStatusCollection> dcsHandle;
     iEvent.getByLabel(dcsTag_, dcsHandle);

     // need the magnetic field
     //
     // if isData then derive bfield using the
     // magnet current from DcsStatus
     // otherwise take it from the IdealMagneticFieldRecord
     if (isData_)
     {
         // scale factor = 3.801/18166.0 which are
         // average values taken over a stable two
         // week period
         float currentToBFieldScaleFactor = 2.09237036221512717e-04;
         float current = -9999/currentToBFieldScaleFactor;
	 if( dcsHandle.isValid() && (*dcsHandle).size() > 0 ) {
	   current = (*dcsHandle)[0].magnetCurrent();
//	   cout << "dcsHandle is valid. The current: " << current << endl;
	 }
	 
         *evt_bField = current*currentToBFieldScaleFactor;
     }
     else
     {
         ESHandle<MagneticField> magneticField;
         iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

         *evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
     }

     std::string branchprefix = aliasprefix_;
     if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

     if( dcsHandle.isValid() && (*dcsHandle).size() > 0 ) {
       *evt_detectorStatus = (*dcsHandle)[0].ready();
       
       iEvent.put(evt_detectorStatus   ,branchprefix+"detectorStatus"  );
     }

     iEvent.put(evt_run              ,branchprefix+"run"             );
     iEvent.put(evt_event            ,branchprefix+"event"           );
     iEvent.put(evt_lumiBlock        ,branchprefix+"lumiBlock"       );
     iEvent.put(evt_bunchCrossing    ,branchprefix+"bunchCrossing"   );
     iEvent.put(evt_orbitNumber      ,branchprefix+"orbitNumber"     );
     iEvent.put(evt_storeNumber      ,branchprefix+"storeNumber"     );
     iEvent.put(evt_experimentType   ,branchprefix+"experimentType"  );
     iEvent.put(evt_bField           ,branchprefix+"bField"          );
     iEvent.put(evt_isRealData       ,branchprefix+"isRealData"      );

}

//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);
