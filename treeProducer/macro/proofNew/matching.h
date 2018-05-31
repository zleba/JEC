//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 27 20:09:14 2018 by ROOT version 6.10/05
// from TTree events/events
// found on file: /nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged/jetsC.root
//////////////////////////////////////////////////////////

#ifndef matching_h
#define matching_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
//#include "KKousour/TopAnalysis/plugins/QCDjet.h"
#include "JEC/treeProducer/interface/QCDjet.h"

#include <vector>
#include <array>

#include "Math/GenVector/PtEtaPhiM4D.h"

#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile2D.h"

#include "../../plugins/JEC.h"


class Luminosity {
    vector<map<int, double>> lumiMap;


    public:

    static map<int,double> GetLumis(TString fName)
    {
        ifstream file(fName);
        assert(file.good());
        
        map<int,double> lumiMap;
        while(1) {
            TString trigName;
            double lum;
            file >> trigName >> lum;
            if(!file.good()) break;
            trigName.ReplaceAll("HLT_PFJet", "");
            trigName = trigName(0, trigName.Length()-3);

            int trigTr = trigName.Atoi();
            if(trigTr != 500) {
                lumiMap[trigTr] += lum;
            }
        }
        file.close();
        return lumiMap;

    }

    void LoadLumis() {
        vector<vector<TString>> files;
        files.push_back({"processedLumis2016B.json.txt"});
        files.push_back({"processedLumis2016C.json.txt"});
        files.push_back({"processedLumis2016D.json.txt"});
        files.push_back({"processedLumis2016E.json.txt"});
        files.push_back({"processedLumis2016Fearly.json.txt", "processedLumis2016Flate.json.txt" }); //Merging F into one period
        files.push_back({"processedLumis2016G.json.txt"});
        files.push_back({"processedLumis2016H.json.txt"});

        lumiMap.resize(files.size() + 1); //Index 0 is dedicated for x-section summed over all periods

        TString path = "/nfs/dust/cms/user/connorpa/SMPJ/effective_luminosity_Run2016BtoH/";

        //Read luminosities to map
        for(unsigned i = 0; i < files.size(); ++i) {
            for(unsigned j = 0; j < files[i].size(); ++j) {
                auto lum = GetLumis(path + files[i][j]);
                for(auto v : lum) {
                    lumiMap[i+1][v.first] += v.second;
                    //cout << "Printout " << v.first <<" "<< v.second << endl;
                }
            }
        }

        //Evaluate the total luminosity
        for(unsigned i = 1; i < lumiMap.size(); ++i)
            for(auto v : lumiMap[i])
                lumiMap[0][v.first] += v.second;

        /*
        //Sum B,C,D and make them identical
        map<int,double> lumiMapSingle;
        for(auto v : lumiMap[0]) {
            double sum = lumiMap[1][v.first] + lumiMap[2][v.first] + lumiMap[3][v.first];
            lumiMap[1][v.first] = lumiMap[2][v.first] = lumiMap[3][v.first] = sum;
        }
        */



        //Assert identical trigger configuration:
        for(unsigned i = 1; i < lumiMap.size(); ++i) {
            for(auto m : lumiMap[0])
                assert(lumiMap[i].count(m.first));
            //assert(lumiMap[0].size() == lumiMap[i].size());
        }


        for(auto v : lumiMap[0]) {
            cout << v.first <<" "<< v.second << endl;
        }


    }

    //runLum, totLum, trigID
    tuple<double,double,int> GetWeightID(int periodId, double pT0)
    {
        //Trigger tresholds from Patrick
        static const map<int,int> trigTrsh = {
            {40 , 74},
            {60 , 84},
            {80 , 114},
            {140 , 196},
            {200 , 245},
            {260 , 330},
            {320 , 395},
            {400 , 468},
            {450 , 507},
        };

        assert(periodId < lumiMap.size());
        int id = trigTrsh.size() - 1;
        for(auto it = trigTrsh.rbegin(); it != trigTrsh.rend(); ++it, --id) 
            if(pT0 >= it->second) {
                return make_tuple(1. / lumiMap[0].at(it->first),  1. / lumiMap[periodId].at(it->first), id);
            }

        return make_tuple(0.,0., -1);

    }

};





class matching : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> runNo = {fReader, "runNo"};
   TTreeReaderValue<Int_t> evtNo = {fReader, "evtNo"};
   TTreeReaderValue<Int_t> lumi = {fReader, "lumi"};
   TTreeReaderValue<Int_t> nvtx = {fReader, "nvtx"};
   TTreeReaderValue<Float_t> pvRho = {fReader, "pvRho"};
   TTreeReaderValue<Float_t> pvz = {fReader, "pvz"};
   TTreeReaderValue<Float_t> pvchi2 = {fReader, "pvchi2"};
   TTreeReaderValue<Float_t> pvndof = {fReader, "pvndof"};
   TTreeReaderValue<Float_t> rho = {fReader, "rho"};
   TTreeReaderValue<Float_t> metEtPF = {fReader, "metEtPF"};
   TTreeReaderValue<Float_t> metSumEtPF = {fReader, "metSumEtPF"};
   TTreeReaderValue<Float_t> metPtPF = {fReader, "metPtPF"};
   TTreeReaderValue<Float_t> metPhiPF = {fReader, "metPhiPF"};
   TTreeReaderValue<Float_t> metEtCHS = {fReader, "metEtCHS"};
   TTreeReaderValue<Float_t> metSumEtCHS = {fReader, "metSumEtCHS"};
   TTreeReaderValue<Float_t> metPtCHS = {fReader, "metPtCHS"};
   TTreeReaderValue<Float_t> metPhiCHS = {fReader, "metPhiCHS"};
   TTreeReaderValue<Float_t> metEtPuppi = {fReader, "metEtPuppi"};
   TTreeReaderValue<Float_t> metSumEtPuppi = {fReader, "metSumEtPuppi"};
   TTreeReaderValue<Float_t> metPtPuppi = {fReader, "metPtPuppi"};
   TTreeReaderValue<Float_t> metPhiPuppi = {fReader, "metPhiPuppi"};

   TTreeReaderValue<std::vector<QCDjet>> chsJets = {fReader, "chsJets"};
   TTreeReaderValue<std::vector<QCDjet>> puppiJets = {fReader, "puppiJets"};

   /*
   TTreeReaderArray<Int_t> chsJets_flavor = {fReader, "chsJets.flavor"};
   TTreeReaderArray<Int_t> chsJets_flavorHadron = {fReader, "chsJets.flavorHadron"};
   TTreeReaderArray<Float_t> chsJets_chf = {fReader, "chsJets.chf"};
   TTreeReaderArray<Float_t> chsJets_nhf = {fReader, "chsJets.nhf"};
   TTreeReaderArray<Float_t> chsJets_phf = {fReader, "chsJets.phf"};
   TTreeReaderArray<Float_t> chsJets_elf = {fReader, "chsJets.elf"};
   TTreeReaderArray<Float_t> chsJets_muf = {fReader, "chsJets.muf"};
   TTreeReaderArray<Int_t> chsJets_chm = {fReader, "chsJets.chm"};
   TTreeReaderArray<Int_t> chsJets_nhm = {fReader, "chsJets.nhm"};
   TTreeReaderArray<Int_t> chsJets_phm = {fReader, "chsJets.phm"};
   TTreeReaderArray<Int_t> chsJets_elm = {fReader, "chsJets.elm"};
   TTreeReaderArray<Int_t> chsJets_mum = {fReader, "chsJets.mum"};
   TTreeReaderArray<Float_t> chsJets_jetJECfact = {fReader, "chsJets.jetJECfact"};
   TTreeReaderArray<Float_t> chsJets_isBtag = {fReader, "chsJets.isBtag"};
   TTreeReaderArray<Float_t> chsJets_btag = {fReader, "chsJets.btag"};
   TTreeReaderArray<Float_t> chsJets_area = {fReader, "chsJets.area"};
   TTreeReaderArray<Float_t> chsJets_unc = {fReader, "chsJets.unc"};
   TTreeReaderArray<Float_t> chsJets_p4_fPt = {fReader, "chsJets.p4.fPt"};
   TTreeReaderArray<Float_t> chsJets_p4_fEta = {fReader, "chsJets.p4.fEta"};
   TTreeReaderArray<Float_t> chsJets_p4_fPhi = {fReader, "chsJets.p4.fPhi"};
   TTreeReaderArray<Float_t> chsJets_p4_fM = {fReader, "chsJets.p4.fM"};
   TTreeReaderArray<Int_t> puppiJets_flavor = {fReader, "puppiJets.flavor"};
   TTreeReaderArray<Int_t> puppiJets_flavorHadron = {fReader, "puppiJets.flavorHadron"};
   TTreeReaderArray<Float_t> puppiJets_chf = {fReader, "puppiJets.chf"};
   TTreeReaderArray<Float_t> puppiJets_nhf = {fReader, "puppiJets.nhf"};
   TTreeReaderArray<Float_t> puppiJets_phf = {fReader, "puppiJets.phf"};
   TTreeReaderArray<Float_t> puppiJets_elf = {fReader, "puppiJets.elf"};
   TTreeReaderArray<Float_t> puppiJets_muf = {fReader, "puppiJets.muf"};
   TTreeReaderArray<Int_t> puppiJets_chm = {fReader, "puppiJets.chm"};
   TTreeReaderArray<Int_t> puppiJets_nhm = {fReader, "puppiJets.nhm"};
   TTreeReaderArray<Int_t> puppiJets_phm = {fReader, "puppiJets.phm"};
   TTreeReaderArray<Int_t> puppiJets_elm = {fReader, "puppiJets.elm"};
   TTreeReaderArray<Int_t> puppiJets_mum = {fReader, "puppiJets.mum"};
   TTreeReaderArray<Float_t> puppiJets_jetJECfact = {fReader, "puppiJets.jetJECfact"};
   TTreeReaderArray<Float_t> puppiJets_isBtag = {fReader, "puppiJets.isBtag"};
   TTreeReaderArray<Float_t> puppiJets_btag = {fReader, "puppiJets.btag"};
   TTreeReaderArray<Float_t> puppiJets_area = {fReader, "puppiJets.area"};
   TTreeReaderArray<Float_t> puppiJets_unc = {fReader, "puppiJets.unc"};
   TTreeReaderArray<Float_t> puppiJets_p4_fPt = {fReader, "puppiJets.p4.fPt"};
   TTreeReaderArray<Float_t> puppiJets_p4_fEta = {fReader, "puppiJets.p4.fEta"};
   TTreeReaderArray<Float_t> puppiJets_p4_fPhi = {fReader, "puppiJets.p4.fPhi"};
   TTreeReaderArray<Float_t> puppiJets_p4_fM = {fReader, "puppiJets.p4.fM"};
   */

   TTreeReaderArray<Float_t> hltJets_fPt = {fReader, "hltJets.fPt"};
   TTreeReaderArray<Float_t> hltJets_fEta = {fReader, "hltJets.fEta"};
   TTreeReaderArray<Float_t> hltJets_fPhi = {fReader, "hltJets.fPhi"};
   TTreeReaderArray<Float_t> hltJets_fM = {fReader, "hltJets.fM"};
   TTreeReaderArray<Float_t> genJets_fPt = {fReader, "genJets.fPt"};
   TTreeReaderArray<Float_t> genJets_fEta = {fReader, "genJets.fEta"};
   TTreeReaderArray<Float_t> genJets_fPhi = {fReader, "genJets.fPhi"};
   TTreeReaderArray<Float_t> genJets_fM = {fReader, "genJets.fM"};
   TTreeReaderValue<vector<bool>> triggerBit = {fReader, "triggerBit"};
   TTreeReaderArray<int> triggerPre = {fReader, "triggerPre"};

    TH1D* hJetPtIncCovTemp,  *hJetPtIncCovTempTot;
    static const int nPer = 8;
    struct Histos {
        TList  *fOutput;
        double w, wTot;
        int fileId;
        array<TH3D*,nPer> hBalEtaPt[5];
        array<TH1D*,nPer> hJetPt;
        array<TH1D*,nPer> hJetPtInc;

        array<TProfile2D*,nPer> hProf[9], hProfBB[9];
                                          
        array<TH3D*,nPer> hRhopuppi;
        array<TH3D*,nPer> hRhochs;
                                          
        array<TH2D*,nPer> hEtaPtCHS;
        array<TH2D*,nPer> hEtaPtPUPPI;
        array<TH2D*,nPer> hEtaPtPUPPIalone;
        void Init(TList *fOutput_);

        void Fill1D(array<TH1D*,nPer> &hist, double x) {
            hist[fileId]->Fill(x, w);
            hist[0]     ->Fill(x, wTot);
        }

        void Fill2D(array<TH2D*,nPer> &hist, double x, double y) {
            hist[fileId]->Fill(x, y, w);
            hist[0]     ->Fill(x, y, wTot);
        }
        void Fill3D(array<TH3D*,nPer> &hist, double x, double y, double z) {
            hist[fileId]->Fill(x, y, z, w);
            hist[0]     ->Fill(x, y, z, wTot);
        }

    } h;

    Luminosity lum;

    JECs jetCorrsCHS, jetCorrsPUPPI;
    TString currFile ="";


    void DoMikkoMatching(int fileId);


   matching(TTree * /*tree*/ =0) { }
   virtual ~matching() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(matching,0);

};

#endif

#ifdef matching_cxx
void matching::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t matching::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef matching_cxx
