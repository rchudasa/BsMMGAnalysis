// -*- C++ -*-
//
// Package:    bs2MuMuGamma/DimuonMassFilter
// Class:      DimuonMassFilter
//
/**\class DimuonMassFilter DimuonMassFilter.cc bs2MuMuGamma/DimuonMassFilter/plugins/DimuonMassFilter.cc

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "TLorentzVector.h"

#define MU_MASS 0.1056583745

#include <vector>
//
// class declaration
//

class DimuonMassFilter : public edm::stream::EDFilter<> {
public:
    explicit DimuonMassFilter(const edm::ParameterSet&);
    ~DimuonMassFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    virtual void beginStream(edm::StreamID) override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

    // input tags
    edm::InputTag SimGenParticleSrc;

    // Tokens
    edm::EDGetTokenT<std::vector<reco::Muon>>         muonToken_;
    

    //
    TLorentzVector mu1_lv,mu2_lv,dimu_lv;
    const bool chargeProd;
    Float_t minptcut;
    Float_t maxAbsEta;
    Float_t minDimuonMass;
    Float_t maxDimuonMass;
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
DimuonMassFilter::DimuonMassFilter(const edm::ParameterSet& iConfig)   :
    chargeProd(iConfig.getUntrackedParameter("ChargeProd", -1)),
    minptcut(iConfig.getUntrackedParameter("MinPt", 3.5)),
    maxAbsEta(iConfig.getUntrackedParameter("MaxAbsEta", 2.5)),
    minDimuonMass(iConfig.getUntrackedParameter("MinDimuonMass", 3.5)),
    maxDimuonMass(iConfig.getUntrackedParameter("MaxDimuonMass", 10.5))
{
    //now do what ever initialization is needed
    muonToken_              = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));

    edm::LogInfo("DimuonMassFilter") << "----------------------------------------------------------------------" << std::endl;
    edm::LogInfo("DimuonMassFilter") << "--- DimuonMassFilter" << std::endl
                                     << " minptcut : " << minptcut<< std::endl
                                     << " maxAbsEta : " << maxAbsEta<< std::endl
                                     << " minDimuonMass :" << minDimuonMass<< std::endl
                                     << " maxDimuonMass :" << maxDimuonMass<< std::endl
                                     << " chargeProd :" << chargeProd<< std::endl;
}


DimuonMassFilter::~DimuonMassFilter()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool DimuonMassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //using namespace edm;
    edm::Handle<std::vector<reco::Muon>> muons;
    iEvent.getByToken(muonToken_, muons);

    bool accepted=false;
    
    //std::cout<<"nMuons : "<<muons->size()<<"\n";
    if(muons->size()<2) return accepted;
    
  for (uint32_t i=0;i<muons->size();i++){
      auto &mu1 = muons->at(i);
      // std::cout<<"mu1 pt eta  : "<<mu1.pt()<<" , "<<mu1.eta()<<" ch : "<<mu1.charge()<<"\n";
      
      if(mu1.pt() < minptcut) continue;
      if(abs(mu1.eta()) > maxAbsEta) continue;
      
     // auto e1= sqrt(mu1.px()*mu1.px() + mu1.py()*mu1.py()+mu1.pz()*mu1.pz() + MU_MASS*MU_MASS);
    //  mu1_lv.SetPxPyPzE( mu1.px(),mu1.py(), mu1.pz(),e1);
      
      for (uint32_t j=i+1;j<muons->size();j++){
        Float_t dimu_mass=0.0;
        auto &mu2 = muons->at(j);
        // std::cout<<"\t mu2 pt eta  : "<<mu2.pt()<<" , "<<mu2.eta()<<" ch : "<<mu2.charge()<<"\n";
        
        if(mu2.charge()*mu2.charge() != chargeProd) continue;
        if(mu2.pt() < minptcut) continue;
        if(abs(mu2.eta()) > maxAbsEta) continue;
        
      //  auto e2= sqrt(mu2.px()*mu2.px() + mu2.py()*mu2.py()+mu2.pz()*mu2.pz() + MU_MASS*MU_MASS);
     //   mu2_lv.SetPxPyPzE( mu2.px(),mu2.py(), mu2.pz(),e2);
        
        dimu_mass=(mu1.p4() + mu2.p4()).mass() ;
        // std::cout<<"\t\t Mass : "<<dimu_mass<<"\n";
        if(dimu_mass < minDimuonMass ) continue;
        if(dimu_mass > maxDimuonMass ) continue;

        accepted= true;
        if(accepted) break;
      }
        if(accepted) break;
    }
   // std::cout<<"Result : "<<accepted<<"\n";
    return accepted;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DimuonMassFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DimuonMassFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DimuonMassFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
DimuonMassFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DimuonMassFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DimuonMassFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DimuonMassFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DimuonMassFilter);
