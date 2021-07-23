// -*- C++ -*-
//
// Package:    bs2MuMuGamma/TriggerFilter
// Class:      TriggerFilter
//
/**\class TriggerFilter TriggerFilter.cc bs2MuMuGamma/TriggerFilter/plugins/TriggerFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Aravind Thachayath Sugunan
//         Created:  Fri, 19 Mar 2021 10:29:17 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// data formats
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"



#include <vector>
//
// class declaration
//

class TriggerFilter : public edm::stream::EDFilter<> {
public:
    explicit TriggerFilter(const edm::ParameterSet&);
    ~TriggerFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginStream(edm::StreamID) override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

  std::string   processName_;
  std::vector<std::string>   sigTriggerNames_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  HLTConfigProvider hltConfig_;

  std::map<std::string,TH1F*> hists_1d_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerFilter::TriggerFilter(const edm::ParameterSet& ps)   
{

  processName_ = ps.getUntrackedParameter<std::string>("processName","HLT");
  std::vector<std::string> trigVec{"HLT_DoubleMu4_3_Bs_v14"};
  sigTriggerNames_ = ps.getUntrackedParameter<std::vector<std::string>>("sigTriggerName",trigVec);
  triggerResultsToken_ = consumes<edm::TriggerResults> (ps.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")));
  verbose_ = ps.getUntrackedParameter<bool>("verbose",false);
    
  // histogram setup
  //edm::Service<TFileService> fs;
  //hists_1d_["h_passSigTrig"] = fs->make<TH1F>("h_passSigTrig" , "; passed ref trigger" , 2 , -0.5 , 1.5 );

 }


TriggerFilter::~TriggerFilter()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool TriggerFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (verbose_) cout << endl;

  // get event products
  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "TriggerFilter::filter: Error in getting TriggerResults product from Event!" << endl;
    return false;
  }
 // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());


  if (verbose_) cout << endl;

  const unsigned int ntrigs(hltConfig_.size());
  unsigned int sigTriggerIndex;
  bool sigAccept=false;
  for(auto& trigName : sigTriggerNames_)
  {
  sigTriggerIndex=hltConfig_.triggerIndex(trigName);
  assert(sigTriggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(trigName));

  // abort on invalid trigger name
   if (sigTriggerIndex>=ntrigs) {
    cout << "TriggerFilter::filterTrigger: path "
	 << trigName << " - not found!" << endl;
    return false;
  }

  if (verbose_) {
    cout << "TriggerFilter::filterTrigger: signal path "
	 << trigName << " [" << sigTriggerIndex << "]" << endl;
  }
  
  // modules on this trigger path
  sigAccept = triggerResultsHandle_->accept(sigTriggerIndex);
  if(verbose_) std::cout<<"trigName : "<<trigName <<" has accept  = :  "<<sigAccept<<"\n"; 
  if(sigAccept) break;
  }
  //if(sigAccept)  hists_1d_["h_passSigTrig"]->Fill(1.0);
  //else		 hists_1d_["h_passSigTrig"]->Fill(0.0);
  
  return sigAccept;
  }

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TriggerFilter::beginStream(edm::StreamID)
{

}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TriggerFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------

void
TriggerFilter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
  for(auto& trigName : sigTriggerNames_)
  {
      unsigned int sigTriggerIndex(hltConfig_.triggerIndex(trigName));
      if (sigTriggerIndex>=n) {
	cout << "TriggerFilter::filter:"
	     << " TriggerName " << trigName
	     << " not available in config!" << endl;
      }
      }
    } // if changed
  } else {
    cout << "TriggerFilter::filter:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }
}


// ------------ method called when ending the processing of a run  ------------
/*
void
TriggerFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TriggerFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TriggerFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TriggerFilter);
