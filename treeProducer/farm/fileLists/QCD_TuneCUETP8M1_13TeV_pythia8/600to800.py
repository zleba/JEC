import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/0048131D-3CB3-E611-813A-001E67DFFB31.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/021372DA-0CB3-E611-80D4-008CFA197DE4.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/06185FAF-00B3-E611-B693-28924A33AFC2.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/06DCD143-03B3-E611-9B4E-0023AEEEB208.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/145C658E-FBB2-E611-87E8-28924A33B9AA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/18836461-3AB3-E611-A3E0-3417EBE47C5E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1A24B5DA-04B3-E611-813B-FA163E05A8CA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/2E1F5AC4-0CB3-E611-945B-28924A3504DA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/324F2D37-06B3-E611-B87F-28924A3504DA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/340AD939-FFB2-E611-B1CA-002590E7DFA2.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3607B5A6-FDB2-E611-BD3B-02163E00E5DC.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3AEB6954-03B3-E611-A16B-C4346BC7EE18.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4410AD6A-FDB2-E611-8067-28924A33AFC2.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/44EB8279-FCB2-E611-ACB5-001E67DFF67C.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/467C5678-3CB3-E611-BB68-009C02AAB554.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/46D7E87E-0CB3-E611-ABBE-0025905A6090.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4E7E1E1A-3CB3-E611-9A2E-0242AC130007.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/50D3AD8F-08B3-E611-A97C-28924A33B062.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/52C4CC80-3CB3-E611-ABFA-001E67792644.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/56B98635-FDB2-E611-8C7B-FA163E8F18B2.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5CA2F6D3-0EB3-E611-99AE-0025905A48D0.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5E782586-3CB3-E611-BE90-001E67F120A8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5E8A7B2D-3CB3-E611-BBDD-FA163EFFA798.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5EFC2431-FAB2-E611-AC88-0090FAA57C60.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/60654951-FAB2-E611-81A4-B8CA3A70A520.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/60F6BF35-FAB2-E611-B3CA-D4AE527EDFE6.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/62466885-00B3-E611-B533-FA163E3CAE59.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/64487D91-08B3-E611-AE88-002590E7DE24.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6E8092CB-3FB3-E611-A42C-008CFA197904.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/725C9A5A-3CB3-E611-B283-02163E0176CD.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/72938A08-03B3-E611-B1FD-FA163E8F0216.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8073365B-FDB2-E611-93B6-90B11C0BD685.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/84BA6AFE-F8B2-E611-8205-002590747DD8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8C11841A-3CB3-E611-8A18-1866DAEA6CC4.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/90AA55C9-00B3-E611-85F8-02163E012F91.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9250CB34-FFB2-E611-9EC1-0023AEEEB208.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/92F98228-3CB3-E611-91EF-001E4F1B8E39.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/962BD9C0-F9B2-E611-B3C5-0CC47AB0B704.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/96D32386-3CB3-E611-8D7D-0025905B85D8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/982B488B-FAB2-E611-9B3E-00304867FCFB.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/98EFD06A-0AB3-E611-9689-008CFA197DDC.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/9CFAE5B1-F9B2-E611-BAEC-28924A33AF26.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/A03781B0-00B3-E611-9CD3-0090FAA57470.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/A67FF63A-08B3-E611-84FA-FA163EC72F40.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/AA3831EB-00B3-E611-9019-008CFA110AA8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/AC93997C-FAB2-E611-AD8E-FA163E08693F.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/AE978146-0BB3-E611-BC06-28924A3504DA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B2911F2F-FFB2-E611-B7E5-00259075D706.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B475DDE2-F7B2-E611-A164-00259029E84C.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B6FB5525-0FB3-E611-A8D1-28924A3504DA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B8D22B32-F9B2-E611-92F0-0025907D2430.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BA1C3E06-F9B2-E611-8D89-02163E013C35.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BA7A910F-07B3-E611-9206-28924A33B062.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BEFEE540-03B3-E611-A78A-0CC47A57D13E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C2201E53-01B3-E611-BFE0-0CC47A57D13E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C2278A36-05B3-E611-B1C8-02163E01316B.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C643C067-FBB2-E611-B27A-002590FD5694.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C8D80E47-FAB2-E611-9CAF-1C6A7A21A11B.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CCC01586-13B3-E611-9107-008CFA1112CC.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CE91E32A-3CB3-E611-AA6C-0CC47AD98C8C.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CEA70CC0-03B3-E611-97A3-002590E7DFA2.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D24B4A47-3CB3-E611-8E48-002590D9D8BE.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D2A09A36-FDB2-E611-A5C0-0CC47A57D13E.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D4B24402-05B3-E611-B393-008CFA111348.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D864F680-3CB3-E611-BD14-00259073E496.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E02BAA00-F9B2-E611-8530-00266CFE6404.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E25E4885-3CB3-E611-9E81-0025905B8568.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E2B33D2F-F8B2-E611-99A1-1866DAEA8218.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E2F7EF20-F7B2-E611-A979-28924A33AFC2.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E4BC6F57-FCB2-E611-A30B-002590FD5A52.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E8C9DB7B-05B3-E611-A9F2-008CFA197D0C.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E8CD1CB1-F9B2-E611-8413-D48564593FA8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/EA2FCA28-FDB2-E611-B98C-002590747DD8.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F0DDA7C6-09B3-E611-8E80-28924A3504DA.root',
'/store/mc/RunIISummer16MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FE4FFA10-10B3-E611-857D-28924A33AFC2.root' ] );


secFiles.extend( [
               ] )