// -*- C++ -*-
//
// Package:    bs2MuMuGamma/GenDecayKineFilter
// Class:      GenDecayKineFilter
//
/**\class GenDecayKineFilter GenDecayKineFilter.cc bs2MuMuGamma/GenDecayKineFilter/plugins/GenDecayKineFilter.cc

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

#include "Pythia8/Pythia.h"

#include <vector>
//
// class declaration
//

class GenDecayKineFilter : public edm::stream::EDFilter<> {
public:
    explicit GenDecayKineFilter(const edm::ParameterSet&);
    ~GenDecayKineFilter();

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
    edm::EDGetTokenT<edm::View<reco::GenParticle>>   SimGenTocken;

    //
    const int fVerbose;
    std::vector<int> dauIDs,antiDauIDs;
    const int particleID;
    int antiParticleID;
    const int motherID;
    const bool chargeconju;
    const int ndaughters;
    std::vector<double> minptcut;
    const double maxptcut;
    std::vector<double> minetacut;
    std::vector<double> maxetacut;
    std::unique_ptr<Pythia8::Pythia> fLookupGen;  // this instance is for accessing particleData information
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
GenDecayKineFilter::GenDecayKineFilter(const edm::ParameterSet& iConfig)   :
    SimGenParticleSrc(iConfig.getParameter<edm::InputTag>("SimGenParticle")),
    SimGenTocken(consumes<edm::View<reco::GenParticle>>(SimGenParticleSrc)),
    fVerbose(iConfig.getUntrackedParameter("verbose", 0)),
    particleID(iConfig.getUntrackedParameter("ParticleID", 0)),
    motherID(iConfig.getUntrackedParameter("MotherID", 0)),
    chargeconju(iConfig.getUntrackedParameter("ChargeConjugation", true)),
    ndaughters(iConfig.getUntrackedParameter("NumberDaughters", 0)),
    maxptcut(iConfig.getUntrackedParameter("MaxPt", 14000.))
{
    //now do what ever initialization is needed

    std::vector<int> defdauID;
    defdauID.push_back(0);
    dauIDs = iConfig.getUntrackedParameter< std::vector<int> >("DaughterIDs", defdauID);
    std::vector<double> defminptcut;
    defminptcut.push_back(0.);
    minptcut = iConfig.getUntrackedParameter< std::vector<double> >("MinPt", defminptcut);
    std::vector<double> defminetacut;
    defminetacut.push_back(-10.);
    minetacut = iConfig.getUntrackedParameter< std::vector<double> >("MinEta", defminetacut);
    std::vector<double> defmaxetacut;
    defmaxetacut.push_back(10.);
    maxetacut = iConfig.getUntrackedParameter< std::vector<double> >("MaxEta", defmaxetacut);

    edm::LogInfo("GenDecayKineFilter") << "----------------------------------------------------------------------" << std::endl;
    edm::LogInfo("GenDecayKineFilter") << "--- GenDecayKineFilter" << std::endl;
    for (unsigned int i = 0; i < dauIDs.size(); ++i) {
        edm::LogInfo("GenDecayKineFilter") << "ID: " << dauIDs[i] << " pT > " << minptcut[i] << " ,  " << minetacut[i]
                                           << "< eta < " << maxetacut[i] << std::endl;
    }
    edm::LogInfo("GenDecayKineFilter") << "maxptcut   = " << maxptcut << std::endl;
    edm::LogInfo("GenDecayKineFilter") << "particleID = " << particleID << std::endl;
    edm::LogInfo("GenDecayKineFilter") << "motherID   = " << motherID << std::endl;

    // create pythia8 instance to access particle data
    edm::LogInfo("GenDecayKineFilter") << "Creating pythia8 instance for particle properties" <<std:: endl;
    if (!fLookupGen.get())
        fLookupGen.reset(new Pythia8::Pythia());

    int antiId;
    if(chargeconju)
    {
        antiParticleID=-particleID;
        if (!(fLookupGen->particleData.isParticle(antiParticleID)))
            antiParticleID = particleID;

        for(size_t i=0; i<dauIDs.size(); i++)
        {
            antiId=-dauIDs[i];
            if (!(fLookupGen->particleData.isParticle(antiId)))
                antiId=dauIDs[i];

            antiDauIDs.push_back(antiId);
        }
    }
}


GenDecayKineFilter::~GenDecayKineFilter()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool GenDecayKineFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle<edm::View<reco::GenParticle>> SimGenCollection;
    iEvent.getByToken(SimGenTocken, SimGenCollection);
    //auto collecA=SimGenCollection.product();

    bool accepted=false;
    std::vector<int> vparticles;
    std::vector<bool> foundDaughter(dauIDs.size(),false);
    std::vector<int> * dauCollection= &dauIDs;
    for(unsigned int i=0; i<SimGenCollection->size(); i++) {

        if((*SimGenCollection)[i].pdgId() == particleID)
        {
            dauCollection = &dauIDs;
        }
        else if(chargeconju and ( (*SimGenCollection)[i].pdgId() == antiParticleID ) )
        {
            dauCollection = &antiDauIDs;
        }
        else {
            continue;
        }
        // -- Check for mother of this particle
        if (0 != motherID) {
            if(((*SimGenCollection)[i].mother(0))->pdgId()!=motherID)
                continue;
        }

        // -- check for daugthers
        int ndau = 0;

        for (unsigned int k = 0; k < foundDaughter.size(); ++k) {
            foundDaughter[k]=false;
        }
        for(unsigned int j=0; j<(*SimGenCollection)[i].numberOfDaughters(); j++)
        {

            auto& daughter = *((*SimGenCollection)[i].daughter(j));

            ++ndau;

            for (unsigned int k = 0; k < dauCollection->size(); ++k) {
                if (daughter.pdgId() != dauCollection->at(k))
                    continue;

                if (daughter.pt() > minptcut[k] && daughter.pt() < maxptcut &&
                        daughter.eta() > minetacut[k] && daughter.eta() < maxetacut[k]) {
                    foundDaughter[k]=true;
                    vparticles.push_back(daughter.pdgId());
                    if (fVerbose > 2) {
                        edm::LogInfo("GenDecayKineFilter")
                                << "  accepted this particle " << daughter.pdgId() << " pT = " << daughter.pt()
                                << " eta = " << daughter.eta() <<std:: endl;
                    }
                    break;
                }
            }
        }


        // -- allow photons
        if (ndau == ndaughters) {
            accepted = true;
            for (unsigned int k = 0; k < foundDaughter.size(); ++k) {
                if(!foundDaughter[k])
                {
                    accepted=false;
                }
            }
            if (accepted and (fVerbose > 0)) {
                edm::LogInfo("GenDecayKineFilter") << "  accepted this decay: ";
                for (unsigned int iv = 0; iv < vparticles.size(); ++iv)
                    edm::LogInfo("GenDecayKineFilter") << vparticles[iv] << " ";
                edm::LogInfo("GenDecayKineFilter") << " from mother = " << motherID <<" -> "<<(*SimGenCollection)[i].pdgId() << std::endl;
            }
        }

        if(accepted) break;
    }

    return accepted;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenDecayKineFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenDecayKineFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenDecayKineFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
GenDecayKineFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenDecayKineFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenDecayKineFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenDecayKineFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenDecayKineFilter);
