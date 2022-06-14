import FWCore.ParameterSet.Config as cms

BmmMuonId = cms.EDProducer(
    "BmmMuonIdProducer",
    muonCollection = cms.InputTag("muons"),
    softMuonMva = cms.FileInPath('BsMMGAnalysis/BsToMuMuGammaNTuplizer/data/muon_mva/Run2018-20210430-2004-Event0.model'),
    isMC = cms.bool(False),
)

BmmMuonIdMc = BmmMuonId.clone( isMC = cms.bool(True) )

def customizeForMuonID(process,pathsToAppend):
    process.BmmMuonId = BmmMuonId
    for pthName in pathsToAppend:
        pth=getattr(process,pthName)
        pth.insert(0,process.BmmMuonId)
