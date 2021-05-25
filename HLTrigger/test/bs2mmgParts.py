process.hltPreDoubleMu43BsToMMG = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtStage2Digis" ),
    offset = cms.uint32( 0 )
)

process.hltDoubleMu4BsToMMGL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    saveTags = cms.bool( True ),
    ChargeOpt = cms.int32( -1 ),
    MaxPtMin = cms.vdouble( 1.0E125 ),
    FastAccept = cms.bool( False ),
    MatchToPreviousCand = cms.bool( True ),
    MaxDr = cms.double( 2.0 ),
    L1CandTag = cms.InputTag( "hltL1fForIterL3L1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0" ),
    inputMuonCollection = cms.InputTag( "hltIterL3Muons" ),
    InputLinks = cms.InputTag( "hltL3MuonsIterL3Links" ),
    PreviousCandIsL2 = cms.bool( True ),
    PreviousCandTag = cms.InputTag( "hltL2fL1sL1DoubleMu0er1p5SQOSdR1p4L1f0L2PreFiltered0" ),
    MaxPtBalance = cms.double( 999999.0 ),
    MaxPtPair = cms.vdouble( 1.0E125 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtMin = cms.vdouble( 3.0 ),
    MaxInvMass = cms.vdouble( 6.0 ),
    MinPtMax = cms.vdouble( 4.0 ),
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    MinN = cms.int32( 1 ),
    MaxDz = cms.double( 9999.0 ),
    MinPtPair = cms.vdouble( 4.9 ),
    CandTag = cms.InputTag( "hltIterL3MuonCandidates" ),
    MinAcop = cms.double( -999.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MinNhits = cms.int32( 0 ),
    NSigmaPt = cms.double( 0.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxEta = cms.double( 2.5 ),
    L1MatchingdR = cms.double( 0.3 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    CutCowboys = cms.bool( False ),
    MinInvMass = cms.vdouble( 4.5 )
)
process.hltDisplacedmumuVtxProducerDoubleMu4BsToMMG = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltIterL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDoubleMu4BsToMMGL3Filtered" ),
    MinPt = cms.double( 0.0 ),
    ChargeOpt = cms.int32( -1 ),
    MaxEta = cms.double( 2.5 ),
    MaxInvMass = cms.double( 999999.0 ),
    MinPtPair = cms.double( 0.0 ),
    matchToPrevious = cms.bool( True ),
    MinInvMass = cms.double( 0.0 )
)
process.hltDisplacedmumuFilterDoubleMu4BsToMMG = cms.EDFilter( "HLTDisplacedmumuFilter",
    saveTags = cms.bool( True ),
    MuonTag = cms.InputTag( "hltIterL3MuonCandidates" ),
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    MinVtxProbability = cms.double( 0.005 ),
    MaxLxySignificance = cms.double( -1.0 ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerDoubleMu4BsToMMG" ),
    FastAccept = cms.bool( True ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinLxySignificance = cms.double( 0.0 )
)

process.HLT_DoubleMu4_3_BsToMMG_v0 = cms.Path( process.HLTBeginSequence + process.hltL1sDoubleMu0er1p5SQOSdRMax1p4IorDoubleMu0er1p4SQOSdRMax1p4 + process.hltPreDoubleMu43BsToMMG + process.hltL1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0 + process.HLTL2muonrecoSequence + cms.ignore(process.hltL2fL1sL1DoubleMu0er1p5SQOSdR1p4L1f0L2PreFiltered0) + process.HLTL3muonrecoSequence + cms.ignore(process.hltL1fForIterL3L1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0) + process.hltDoubleMu4BsToMMGL3Filtered + process.hltDisplacedmumuVtxProducerDoubleMu4BsToMMG + process.hltDisplacedmumuFilterDoubleMu4BsToMMG + process.HLTEndSequence )


