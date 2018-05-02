#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
//#include "Math/SpecFuncMathMore.h"
#include "TMath.h"
//#include "TVectorD.h"
//#include "TMatrixDSym.h"
//#include "TMatrixDSymEigen.h"

//#include "fastjet/contrib/Njettiness.hh"
//#include "fastjet/tools/MassDropTagger.hh"
//#include "fastjet/contrib/SoftDrop.hh"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "JEC/treeProducer/plugins/treeProducer.h"


#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"


#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


using namespace std;
using namespace reco;
//using namespace fastjet;

treeProducer::treeProducer(edm::ParameterSet const& cfg) : p(cfg, consumesCollector())
{ }



/*
vector<float> getArray(vector<QCDjet> &jets)
{

}
*/

vector<QCDjet> FillJets (JECs &jetEcorrs, edm::Handle<pat::JetCollection> &jets, double rho)
{
  vector<QCDjet> jetVec;
  for(pat::JetCollection::const_iterator ijet =jets->begin();ijet != jets->end(); ++ijet) {
      //if(ijet->pt() < 20) continue;
      //if (isGoodJet(*ijet)) {
      QCDjet jetNow;
      jetNow = GetJet(*ijet);

      vector<string> dumy;
      double L2L3res, Unc;
      //cout << "PT before " << ijet->pt() << endl;
      //cout << "Rho is " << *rho <<" "<<  pvRho_  <<  endl;
      jetNow.jetJECtot = 
           jetEcorrs.JEC_CHScorrections(ijet->correctedP4("Uncorrected").pt(), ijet->eta(), ijet->jetArea(),  rho, dumy, L2L3res, Unc);
      jetNow.jetJECl2l3Res = L2L3res;
      //cout << "PT after " << jet.Pt() << endl;
      //jetNow.p4 =  ROOT::Math::PtEtaPhiM4D<float>(newPt, ijet->eta(), ijet->phi(), ijet->mass()); 

      jetNow.btag = ijet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"/*  p.srcBtag_.c_str()*/);

      //mPFUncCHS.setJetEta(ijet->eta());
      //mPFUncCHS.setJetPt(jet.Pt()); // here you must use the CORRECTED jet pt
      jetNow.unc = Unc;
      //cout << "Jet unc is " << jet.Pt()<<" "<< Unc << endl;
      //jetNow.jetJECfact = CorFactorL2L3res;

      //push only jets above 10 GeV
      if(jetNow.p4.Pt() < 10) continue;

      jetVec.push_back(jetNow);
  }

  sort(jetVec.begin(), jetVec.end(), [](QCDjet &a, QCDjet &b){return  a.p4.Pt() > b.p4.Pt();});

  return jetVec;

}





//////////////////////////////////////////////////////////////////////////////////////////
void treeProducer::beginJob() 
{
  cout << "Starting treeProducer " << __LINE__ << endl;
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetCanExtend(TH1::kAllAxes);
  for(unsigned i=0;i<p.triggerNames_.size();i++) {
    triggerNamesHisto_->Fill(p.triggerNames_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetCanExtend(TH1::kAllAxes);

  cutFlowHisto_ = fs_->make<TH1F>("CutFlow","CutFlow",1,0,1);
  cutFlowHisto_->SetCanExtend(TH1::kAllAxes);


  //Init of JEC

  char period = 'P';
  if(!p.isMC_) {
      const string sPatern = "data/Run2016";
      cout << "RADEK " << p.curFile_ << endl;
      size_t place = p.curFile_.find(sPatern);
      assert(place != string::npos);
      cout << p.curFile_ << endl;
      period = p.curFile_[place + sPatern.size()];
      cout << "period is " << period << endl;
  }



  string jecTag = "Summer16_07Aug2017";
  int version = 7;
  //string jecTag = "Spring16_25ns";
  //int version = 6;


  string jetType = "AK4PFchs";
  if(p.jetType_.find("Puppi") != string::npos) {
      jetType = "AK4PFPuppi";
      //cout << "Puppi found" << endl;
      //exit(0);
  }

  vector<string> dumy;

  cout << "Jet type is :"  << p.isMC_ <<" "<< jetType << endl;
  jetEcorrsCHS.Init(p.isMC_, jecTag, period, version, "AK4PFchs", "", dumy);
  jetEcorrsPUPPI.Init(p.isMC_, jecTag, period, version, "AK4PFPuppi", "", dumy);
  //exit(0);
 
  //--- book the tree ----------------------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"runNo/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evtNo/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nvtx/I");
  outTree_->Branch("pvRho"                ,&pvRho_             ,"pvRho/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz/F");
  outTree_->Branch("pvchi2"               ,&pvchi2_            ,"pvchi2/F");
  outTree_->Branch("pvndof"               ,&pvndof_            ,"pvndof/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho/F");

  //outTree_->Branch("ht"                   ,&ht_                ,"ht/F");
  //
  outTree_->Branch("metEtPF"                ,&metEtPF_             ,"metEtPF/F");
  outTree_->Branch("metSumEtPF"             ,&metSumEtPF_          ,"metSumEtPF/F");
  outTree_->Branch("metPtPF"               ,&metPtPF_             ,"metPtPF/F");
  outTree_->Branch("metPhiPF"              ,&metPhiPF_            ,"metPhiPF/F");

  outTree_->Branch("metEtCHS"            ,&metEtCHS_         ,"metEtCHS/F");
  outTree_->Branch("metSumEtCHS"         ,&metSumEtCHS_      ,"metSumEtCHS/F");
  outTree_->Branch("metPtCHS"           ,&metPtCHS_         ,"metPtCHS/F");
  outTree_->Branch("metPhiCHS"          ,&metPhiCHS_        ,"metPhiCHS/F");

  outTree_->Branch("metEtPuppi"           ,&metEtPuppi_        ,"metEtPuppi/F");
  outTree_->Branch("metSumEtPuppi"        ,&metSumEtPuppi_     ,"metSumEtPuppi/F");
  outTree_->Branch("metPtPuppi"          ,&metPtPuppi_        ,"metPtPuppi/F");
  outTree_->Branch("metPhiPuppi"         ,&metPhiPuppi_       ,"metPhiPuppi/F");

  outTree_->Branch("chsJets"           ,&chsJets_);
  outTree_->Branch("puppiJets"         ,&puppiJets_);


  //  outTree_->Branch("metPuppimass_"        ,&metPuppimass_      ,"metPuppimass_/F");
  //------------------------------------------------------------------

  HLTjets_        = new std::vector< ROOT::Math::PtEtaPhiM4D<float> >;

  chsJets_         = new std::vector<QCDjet>;
  puppiJets_       = new std::vector<QCDjet>;
  genJets_         = new std::vector<ROOT::Math::PtEtaPhiM4D<float>>;

  /*
  */
  //------------------------------------------------------------------
  outTree_->Branch("hltJets"           ,  &HLTjets_);
  outTree_->Branch("genJets"           ,  &genJets_);


  //------------------------------------------------------------------
  triggerBit_ = new std::vector<bool>;
  triggerPre_ = new std::vector<int>;
  outTree_->Branch("triggerBit"           ,"vector<bool>"      ,&triggerBit_);
  outTree_->Branch("triggerPre"           ,"vector<int>"       ,&triggerPre_);

  //outTree_->Branch("nTriggerObjects", &nTriggerObjects_, "nTriggerObjects/I");

  cout<<"RADEK Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void treeProducer::endJob() 
{  
}
//////////////////////////////////////////////////////////////////////////////////////////
void treeProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
  if (p.isMC_ && p.debug_) {
    iRun.getByToken(p.runInfoToken,runInfo);
    for(vector<LHERunInfoProduct::Header>::const_iterator it = runInfo->headers_begin();it != runInfo->headers_end(); it++) {
      cout<<it->tag()<<endl;
      vector<string> lines = it->lines();
      for(unsigned int iLine = 0; iLine < lines.size(); iLine++) {
	cout<< lines.at(iLine);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void treeProducer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
bool treeProducer::isGoodJet(const pat::Jet &jet)
{
  bool res  = true; // by default is good, unless fails a cut bellow
/*  float chf = jet.chargedHadronEnergyFraction();
  float nhf = jet.neutralHadronEnergyFraction();
  float phf = jet.photonEnergyFraction();
  float muf = jet.muonEnergyFraction();
  float elf = jet.electronEnergyFraction();
  int chm   = jet.chargedHadronMultiplicity();
  int npr   = jet.neutralMultiplicity()+jet.chargedMultiplicity();
  float eta = fabs(jet.eta());
  float pt  = jet.pt();
  bool idL  = (npr>1 && phf<0.99 && nhf<0.99);
  //bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
  bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || eta>2.4));
  //if (!idT) res = false;
  if (pt < ptMin_) res = false;
  if (eta > etaMax_) res = false;
  cout << "Jet Parameters "<< idT << " "<< pt <<" "<< eta << endl;
  //if (jet.userFloat("ak8PFJetsCHSSoftDropMass") < massMin_) res = false;
  //if ((jet.subjets("SoftDrop")).size() < 2) res = false;*/
  return res;
}
//////////////////////////////////////////////////////////////////////////////////////////
void treeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{

  initialize();

  //  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;

  //    if(p.isPrint_) cout<<"**** EVENT ****"<<endl;

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects1;

  iEvent.getByToken(p.genjetsToken, genjets);
  iEvent.getByToken(p.jetsCHSToken,jetsCHS);
  iEvent.getByToken(p.jetsPUPPIToken,jetsPUPPI);
  iEvent.getByToken(p.candsToken,cands);
  iEvent.getByToken(p.rhoToken,rho);
  iEvent.getByToken(p.recVtxsToken,recVtxs);  
  iEvent.getByToken(p.triggerResultsToken,triggerResults);  
  iEvent.getByToken(p.triggerPrescalesToken,triggerPrescales); 
  iEvent.getByToken(p.triggerObjectsToken, triggerObjects1);
  iEvent.getByToken(p.met1Token,met1);
  iEvent.getByToken(p.met2Token,met2);
  iEvent.getByToken(p.met3Token,met3);
  //iEvent.getByToken(p.metCHSToken,metCHS);

  triggerBit_->clear();
  triggerPre_->clear();

  //-------------- Trigger Info -----------------------------------
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);  

  //  assert(triggerObjects1);
  // cout << "TrigArr val " <<  triggerObjects << endl;
  //cout << "Trig size " << triggerObjects->size() << endl;
  //for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      //cout << "hi Radek" << endl;
      //obj.unpackPathNames(names);

      //std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      //std::vector<std::string> pathNamesLast = obj.pathNames(true);

      //cout << "RADEK " << obj.pt() <<" "<< pathNamesAll.size() << " "<< pathNamesLast.size() << endl;
      //TLorentzVector P4;                                                                                                                                         
      //P4.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());                                                                                                      
      //LorentzVector qcdhltobj(P4.Px(),P4.Py(),P4.Pz(),P4.E());             
  //}


  set<string> trigNames;
  triggerPassHisto_->Fill("totalEvents",1);
  bool passTrigger(false);
  for(unsigned int k=0; k<p.triggerNames_.size(); ++k) {
      bool bit(false);
      int pre(1);
      for(unsigned int itrig=0; itrig<triggerResults->size(); ++itrig) {
          //cout<<"Trigger name of index "<<itrig<<" "<<string(names.triggerName(itrig))<<endl;
          string trigger_name = string(names.triggerName(itrig));
          //--- erase the last character, i.e. the version number----
          trigger_name.pop_back();
          if (trigger_name == p.triggerNames_[k]) {
              bit = triggerResults->accept(itrig);
              if(bit) trigNames.insert(trigger_name);
              pre = triggerPrescales->getPrescaleForIndex(itrig);
              if (bit) {
                  //cout << "Passed " << p.triggerNames_[k].c_str() << endl;
                  triggerPassHisto_->Fill(p.triggerNames_[k].c_str(),1);
              } 
          }
      }
      //--- if at least one monitored trigger has fired passTrigger becomes true
      passTrigger += bit;
      triggerBit_->push_back(bit); 
      triggerPre_->push_back(pre);

  }   

  //skip events without fired trigger
  //cout << "passTriggers " << passTrigger <<" " << nTrue << endl;
  if(!p.isMC_ && !passTrigger) return;

  //exit(0);
  /*
  set<string> trigNames;
  //for(unsigned int k=0; k<p.triggerNames_.size(); ++k) {
  for(unsigned int itrig=0; itrig<triggerResults->size(); ++itrig) {
      string trigger_name = string(names.triggerName(itrig));
      trigger_name.pop_back();

      bool isInOurTriggs = false;
      for(unsigned int k=0; k<p.triggerNames_.size(); ++k) 
          if (trigger_name == p.triggerNames_[k]) {isInOurTriggs = true; break; }

      if(!isInOurTriggs) continue;

      trigNames.insert(trigger_name);
      cout << "Triggers which fired " <<itrig <<" "<< trigger_name << endl;
  }


  cout << "My nice triggers start" << endl;
  for(auto a : trigNames)
      cout << a << endl;
  cout << "My nice triggers end" << endl;
  */


  nTriggerObjects_ = 0;
  //cout<<"let find HLTobject"<<endl;

  //  cout << "Starting HLT loop" << endl;
  vector<ROOT::Math::PtEtaPhiM4D<float> > hltVecs;
  for(pat::TriggerObjectStandAlone obj: *triggerObjects1){
      obj.unpackPathNames(names);
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);                                                                            
      std::vector<std::string> pathNamesLast = obj.pathNames(true);

      if(pathNamesAll.size() == 1 && pathNamesLast.size() == 1) {
          string nTemp = pathNamesAll[0];
          nTemp.pop_back();

          if(trigNames.count(nTemp) > 0) {
              //cout <<"pT " <<setprecision(7) <<  obj.pt()<<"\t" << obj.eta() <<" "<< obj.phi() <<" "<< nTemp << endl;
              ROOT::Math::PtEtaPhiM4D<float> P4(obj.pt(), obj.eta(), obj.phi(), obj.mass());
              //P4.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());
              bool isIn = false;
              for(const auto &v : hltVecs)
                  if(v == P4) {isIn = true; break;}
              if(!isIn){
                  hltVecs.push_back(P4);
              }
          }
      }
      /*
      bool isOurTrigger = false;
      for(unsigned j=0; j<pathNamesAll.size();j++) {
          string nTemp = pathNamesAll[j];
          nTemp.pop_back();
          if(trigNames.count(nTemp) > 0) {
              isOurTrigger = true;
              break;
          }
      }
      if(!isOurTrigger) continue;
      //cout << "Now is used " << isOurTrigger << endl;
      cout <<"pT " <<setprecision(7) <<  obj.pt()<<"\t" << obj.eta() <<" "<< obj.phi() <<" "<< pathNamesAll.size()<<"\t"<<pathNamesLast.size()<<endl;    
      for(unsigned j=0; j<pathNamesAll.size();j++) cout<<"all\t"<<pathNamesAll[j]<<endl;
      //for(unsigned j=0; j<pathNamesLast.size();j++) cout<<"last\t"<<pathNamesLast[j]<<endl;

      */
      //      ++nTriggerObjects_;
  }
  std::sort(hltVecs.begin(), hltVecs.end(), [](const ROOT::Math::PtEtaPhiM4D<float> &v1, const ROOT::Math::PtEtaPhiM4D<float> &v2) { return v1.Pt() > v2.Pt(); });
  for(const auto &v : hltVecs) {
      HLTjets_->push_back((ROOT::Math::PtEtaPhiM4D<float>(v.Pt(), v.Eta(), v.Phi(), v.M())  ));
  }

  /*
  for(const auto &v : hltVecs)
    cout <<"HLT " <<  v.Pt() << " "<< v.Eta()  << " "<< v.Phi() << endl;
  */

  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  if (cut_vtx) {
    pvRho_    = (*recVtxs)[0].position().Rho();
    pvz_    = (*recVtxs)[0].z();
    pvndof_ = (*recVtxs)[0].ndof();
    pvchi2_ = (*recVtxs)[0].chi2();
  }// if vtx
  //----- PF jets ------------------------------


  //edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParCollCHS;
  //iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",PFJetCorParCollCHS);
  //JetCorrectorParameters const& PFJetCorParCHS = (*PFJetCorParCollCHS)["Uncertainty"];

  //mPFUncCHS = new JetCorrectionUncertainty(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");
  //JetCorrectionUncertainty mPFUncCHS(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");





  vector<QCDjet> jetVecCHS = FillJets(jetEcorrsCHS, jetsCHS, *rho);
  chsJets_ = &jetVecCHS;
  vector<QCDjet> jetVecPUPPI = FillJets(jetEcorrsPUPPI, jetsPUPPI, *rho);
  puppiJets_ = &jetVecPUPPI;



  /*
  for(pat::JetCollection::const_iterator ijet =jetsCHS->begin();ijet != jetsCHS->end(); ++ijet) {
      auto mySets   = ijet->availableJECSets();
      auto myLevels = ijet->availableJECLevels();
      vector<string> dumy;
      double L2L3res, Unc;
      double corr = jetEcorrsCHS.JEC_CHScorrections(ijet->correctedP4("Uncorrected").pt(), ijet->eta(), ijet->jetArea(),  *rho, dumy, L2L3res, Unc);
      cout << "Radek " <<  ijet->pt() / ijet->correctedP4("Uncorrected").pt()  << " "<< corr << endl;
      //cout << ijet->correctedP4("Uncorrected").pt()
  }
  */
















  //vector<LorentzVector> vP4Gen;


  if(p.isMC_) {
      for(GenJetCollection::const_iterator igen = genjets->begin(); igen != genjets->end(); ++igen) {
          genJets_->push_back(ROOT::Math::PtEtaPhiM4D<float>(igen->pt(),igen->eta(),igen->phi(),igen->mass()) );
      }
  }




  rho_    = *rho;
  nVtx_   = recVtxs->size();
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();

  metEtPF_ = (*met1)[0].et();
  metSumEtPF_ = (*met1)[0].sumEt();
  metPtPF_ = (*met1)[0].pt();
  metPhiPF_ = (*met1)[0].phi();

  metEtCHS_ = (*met2)[0].et();
  metSumEtCHS_ = (*met2)[0].sumEt();
  metPtCHS_ = (*met2)[0].pt();
  metPhiCHS_ = (*met2)[0].phi();
  //  metNoHFmass_ = (*met2)[0].mass();

  metEtPuppi_ = (*met3)[0].et();
  metSumEtPuppi_ =(*met3)[0].sumEt();
  metPtPuppi_ = (*met3)[0].pt();
  metPhiPuppi_ = (*met3)[0].phi();


  //cout <<"Met types " << metEtPF_ << " "<<  metEtCHS_ <<" "<< metEtPuppi_ << endl;

  //cout << metCHS->size()<< "  " << metEt_ << " "<< (*metCHS)[0].et() <<  endl;

  //metPuppimass_ = (*met3)[0].mass();

  //cout << "I am filling " << endl;
  outTree_->Fill();



}



//////////////////////////////////////////////////////////////////////////////////////////
void treeProducer::initialize()
{
  run_            = -1;
  evt_            = -1;
  lumi_           = -1;
  nVtx_           = -1;
  nJets_          = -1;
  rho_            = -1;
  metEtPF_        = -1;
  metEtCHS_      = -1;
  metEtPuppi_     = -1;
  metSumEtPF_       = -1;
  metSumEtCHS_   = -1;
  metSumEtPuppi_  = -1;
  metPtPF_          = -1;
  metPhiPF_         = -1;
  metPtCHS_      = -1;
  metPhiCHS_     = -1;
  metPtPuppi_     = -1;
  metPhiPuppi_    = -1;
  //metPuppimass_   = -1;
  pvRho_          = -999;
  pvz_            = -999;
  pvndof_         = -999;
  pvchi2_         = -999;
  genJets_        ->clear();
  HLTjets_        ->clear();
  //qcdJet_->clear();

  nTriggerObjects_ = -1;
}
//////////////////////////////////////////////////////////////////////////////////////////
treeProducer::~treeProducer() 
{
}

DEFINE_FWK_MODULE(treeProducer);
















