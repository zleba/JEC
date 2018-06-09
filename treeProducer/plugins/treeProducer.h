#ifndef treeProducer_h
#define treeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/Math/interface/deltaPhi.h"
//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"//add
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"//add

//#include "JEC.h"
#include "JEC/treeProducer/interface/QCDjet.h"


#include <sstream>
/*
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
*/

using namespace reco;

inline QCDjet GetJet(const pat::Jet &ijet) {
    QCDjet jjet;
    jjet.flavor        =ijet.partonFlavour();
    jjet.flavorHadron  =ijet.hadronFlavour();
    jjet.chf           =ijet.chargedHadronEnergyFraction();
    jjet.nhf           =ijet.neutralHadronEnergyFraction(); //NHF
    jjet.phf           =ijet.photonEnergyFraction();
    jjet.elf           =ijet.electronEnergyFraction();
    jjet.muf           =ijet.muonEnergyFraction(); //MUF
    jjet.chm           =ijet.chargedHadronMultiplicity();
    jjet.nhm           =ijet.neutralHadronMultiplicity();
    jjet.phm           =ijet.photonMultiplicity();
    jjet.elm           =ijet.electronMultiplicity();
    jjet.mum           =ijet.muonMultiplicity();
    jjet.area          =ijet.jetArea();
    //jjet.p4            =ROOT::Math::PtEtaPhiM4D<float>( ijet.correctedP4("Uncorrected").pt(), ijet.eta(), ijet.phi(), ijet.correctedP4("Uncorrected").mass() );
    jjet.p4            = ijet.correctedP4("Uncorrected");
    return jjet;
}




struct Parameters {
    edm::EDGetTokenT<pat::JetCollection> jetsCHSToken, jetsPUPPIToken;
    edm::EDGetTokenT<GenJetCollection> genjetsToken;
    //  muonsToken         
    //electronsToken       
    edm::EDGetTokenT<pat::METCollection> met1Token;
    edm::EDGetTokenT<pat::METCollection> met2Token;
    edm::EDGetTokenT<pat::METCollection> met3Token;
    //edm::EDGetTokenT<pat::METCollection> metCHSToken;

    edm::EDGetTokenT<pat::PackedCandidateCollection> candsToken;
    edm::EDGetTokenT<double> rhoToken;
    edm::EDGetTokenT<reco::VertexCollection> recVtxsToken;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pupInfoToken;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken;
    edm::EDGetTokenT<LHEEventProduct> lheEvtInfoToken;
    edm::EDGetTokenT<LHERunInfoProduct> runInfoToken;
    std::vector<std::string> triggerNames_;
    std::string srcBtag_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
  

    HLTConfigProvider hltConfig_;

    bool isMC_ ;
    bool isPrint_ ;
    bool saveWeights_;
    bool debug_ ;
    string jetType_;
    string curFile_;

    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

    Parameters(edm::ParameterSet const& cfg,  edm::ConsumesCollector && iC) {
        jetsCHSToken             = iC.consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jetsCHS"));
        jetsPUPPIToken           = iC.consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jetsPUPPI"));

        stringstream temp; temp << cfg.getParameter<edm::InputTag>("jetsCHS");
        jetType_ = temp.str();

        genjetsToken          = iC.consumes<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")));
        //  muonsToken            = iC.consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"));
        //electronsToken        = iC.consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("electrons"));
        met1Token              = iC.consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met1"));
        met2Token              = iC.consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met2"));
        met3Token              = iC.consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met3"));

        //metCHSToken            = iC.consumes<pat::METCollection>(string("slimmedMETsCHS"));
        candsToken            = iC.consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("candidates"));
        rhoToken              = iC.consumes<double>(cfg.getParameter<edm::InputTag>("rho"));
        recVtxsToken          = iC.consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"));
        triggerResultsToken   = iC.consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"));
        triggerPrescalesToken = iC.consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("triggerPrescales"));
        pupInfoToken          = iC.consumes<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
        genEvtInfoToken       = iC.consumes<GenEventInfoProduct>(edm::InputTag("generator"));
        genParticlesToken     = iC.consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
        lheEvtInfoToken       = iC.consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
        runInfoToken          = iC.consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));
        srcBtag_              = cfg.getParameter<std::string>("btagger");
        //  xmlFile_              = cfg.getParameter<std::string>("xmlFile");
        triggerNames_         = cfg.getParameter<std::vector<std::string> >("triggerNames");

        triggerObjectsToken  = iC.consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag>("triggerObjects"));

        //massMin_              = cfg.getParameter<double>("massMin");
        //minMuPt_              = cfg.getParameter<double>("minMuPt");
        //minElPt_              = cfg.getParameter<double>("minElPt");
        isMC_                 = cfg.getUntrackedParameter<bool>("isMC",false);
        isPrint_              = cfg.getUntrackedParameter<bool>("isPrint",false);
        saveWeights_          = cfg.getUntrackedParameter<bool>("saveWeights",true);
        debug_                = cfg.getUntrackedParameter<bool>("debug",false);

        //cout <<"RADEK listFile " <<  cfg.getParameter<edm::InputTag>("listFile") << endl;

        //curFile_      =           cfg.getParameter<string>("fileNames");

        jetFlavourInfosToken_ = iC.consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));
    }
};






class treeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit treeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~treeProducer();

  private:  
//    virtual bool isGoodMuon(const pat::Muon &mu,edm::Handle<pat::PackedCandidateCollection> pfcands);
//    virtual bool isGoodElectron(const pat::Electron &el,const reco::Vertex &vtx,edm::Handle<pat::PackedCandidateCollection> pfcands);
    virtual bool isGoodJet(const pat::Jet &jet);
//    float getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,const reco::Candidate *cand);
    void initialize();

    //---- configurable parameters --------  
    Parameters p;

    //    edm::EDGetTokenT<pat::MuonCollection> muonsToken;
    //edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
 
    //std::string srcBtag_;//,xmlFile_;

    //---------------------------------------------------
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    TH1F *cutFlowHisto_;
    //---- TRIGGER -------------------------  
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int   run_,evt_,nVtx_,lumi_,nJets_,nBJets_,nLeptons_,nGenJets_,nTriggerObjects_;
    float rho_,met_,metSig_,ht_,mva_,pvRho_,pvz_,pvndof_,pvchi2_,mvaGen_,metGenSig_;
    float metEtCHS_,metSigEt_,metSumEtCHS_,metEtPF_,metSigEtPF_,metSumEtPF_,metEtPuppi_,metSigEtPuppi_,metSumEtPuppi_;
    float metPtCHS_,metPhiCHS_,metPtPF_,metPhiPF_,metPtPuppi_,metPhiPuppi_;
    vector<QCDjet> *chsJets_, *puppiJets_;
    vector< ROOT::Math::PtEtaPhiM4D<float> > *HLTjets_, *genJets_;

    std::vector<bool> *triggerBit_;
    std::vector<int>  *triggerPre_;
    //---- top variables --------------
    float dRJJ_,dPhiJJ_,mJJ_,yJJ_,ptJJ_;
    float dRGenJJ_,dPhiGenJJ_,mGenJJ_,yGenJJ_,ptGenJJ_,metGen_;
    float dPhiLJ_;
    //---- jet variables --------------
    std::vector<int>   *flavorSub0_,*flavorSub1_,*flavorHadronSub0_,*flavorHadronSub1_;

    //----- HLT jet variables
    std::vector<float> *HLTpt_,*HLTeta_,*HLTphi_,*HLTmass_;

    //---- lepton variables -----------
    std::vector<int>   *lId_;
    std::vector<float> *lPt_,*lEta_,*lPhi_,*lE_,*lIso_;

    std::vector<int>   *lGenId_;
    std::vector<float> *lGenPt_,*lGenEta_,*lGenPhi_,*lGenE_;
    
    std::vector<float> *GenSubJet1Pt_,*GenSubJet2Pt_,*GenSubJet1Eta_,*GenSubJet2Eta_,*GenSubJet1Phi_,*GenSubJet2Phi_,*GenSubJet1Mass_,*GenSubJet2Mass_,*GenSubJetsDeltaR_,*GenSubJetsMu_;
    //---- MVA discriminator ----------
    //BoostedDiscriminatorMVA *discr_;
    //BoostedDiscriminatorMVA *discrGen_;
    //---- MC variables ---------------
    int npu_,decay_;
    float genEvtWeight_,lheOriginalXWGTUP_;
    std::vector<float> *scaleWeights_;
    std::vector<float> *pdfWeights_;
    float ptTopParton_[2],yTopParton_[2],mTTbarParton_,yTTbarParton_,ptTTbarParton_;
    std::vector<int>   *partonId_,*partonSt_,*partonMatchIdx_;
    std::vector<float> *partonPt_,*partonEta_,*partonPhi_,*partonE_,*partonMatchDR_;

    std::vector<int>   *WBosonId_,*WBosonSt_;
    std::vector<float> *WBosonPt_,*WBosonEta_,*WBosonPhi_,*WBosonE_;

    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

    //gen jets
    std::vector<float> *GenSoftDropMass_,*GenSoftDropTau32_,*GenSoftDropTau31_;
    std::vector<float> *GenJetpt_;
    std::vector<float> *GenJetphi_;
    std::vector<float> *GenJeteta_;
    std::vector<float> *GenJetenergy_;
    std::vector<float> *GenJetmass_;
    std::vector<bool> *isBJetGen_;
    std::vector<bool> *isWJetGen_;

    edm::Handle<pat::JetCollection> jetsCHS, jetsPUPPI;
    edm::Handle<GenJetCollection> genjets;
    //edm::Handle<pat::MuonCollection> muons;
    //edm::Handle<pat::ElectronCollection> electrons;
    edm::Handle<pat::METCollection> met1;
    edm::Handle<pat::METCollection> met2;
    edm::Handle<pat::METCollection> met3;
    //edm::Handle<pat::METCollection> metCHS;
    edm::Handle<pat::PackedCandidateCollection> cands;
    edm::Handle<double> rho;
    edm::Handle<reco::VertexCollection> recVtxs;
    edm::Handle<edm::TriggerResults> triggerResults;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::TriggerObjectStandAlone> triggerObjects;
    edm::Handle<edm::View<PileupSummaryInfo> > pupInfo;
    edm::Handle<edm::View<reco::GenParticle> > genParticles;
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    edm::Handle<LHEEventProduct> lheEvtInfo;
    edm::Handle<LHERunInfoProduct> runInfo;

    //JetCorrectionUncertainty *mPFUncCHS;
    //JECs jetEcorrsCHS, jetEcorrsPUPPI;

    //fastjet                                                                                                                                                                   
    //fastjet::Filter* fTrimmer1;
    //fastjet::JetDefinition*       fAKJetDef;
    //fastjet::ActiveAreaSpec*      fActiveArea;
    //fastjet::AreaDefinition*      fAreaDefinition;
    //fastjet::ClusterSequenceArea* fClustering;
    //fastjet::contrib::SoftDrop* sd;
};





#endif
