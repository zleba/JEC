#define matching_cxx
// The class definition in matching.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("matching.C")
// root> T->Process("matching.C","some options")
// root> T->Process("matching.C+")
//


#include "matching.h"
#include <TH2.h>
#include <TStyle.h>
#include "TRandom.h"


//const vector<double> etaBins2  = { 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

const vector<double> etaBinsPos =  {0., 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.538, 4.716, 4.889, 5.191};


const vector<double> Ptbinning = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,  2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717,    7000};






//double dist2(double eta1, double phi1, double eta2, double phi2)
double dist2(const ROOT::Math::PtEtaPhiM4D<float> &j1, const ROOT::Math::PtEtaPhiM4D<float> &j2)
{
    double dPhi = abs(j2.Phi() - j1.Phi());
    dPhi = min(dPhi, 2*M_PI - dPhi);
    double dEta = abs(j2.Eta() - j1.Eta());

    return (dEta*dEta + dPhi*dPhi);
}

const bool withoutRes = false; //without residual corrections?

const string jecTagCHS = "Summer16_07Aug2017";
const int versionCHS = 11;

const string jecTagPUPPI = "Summer16_07Aug2017";
const int versionPUPPI = 11;

//const string jecTagPUPPI = "Spring16_23Sep2016";
//const int versionPUPPI = 2;


bool isLeptonVetoJet(const QCDjet &jet) {

    auto NHF = jet.nhf; //OK
    auto NEMF = jet.phf;//OK
    auto NumConst = jet.chm + jet.nhm + jet.phm + jet.elm + jet.mum; //OK

    auto MUF = jet.muf;//OK
    auto eta = jet.p4.Eta();//OK
    auto CHF = jet.chf;//OK
    auto CHM = jet.chm;//OK
    auto CEMF= jet.elf + jet.muf;//OK

    bool tightLepVetoJetID = true;

    if(abs(eta)<=2.7)
       tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4);
    
    return tightLepVetoJetID;
}



//#ifdef __MAKECINT__
//#pragma link C++ class vector<Jet>+;
//#endif

//gInterpreter->GenerateDictionary("Jet","match.h");
//#include "match.h"


//#pragma link C++ class vector<Jet> +;

struct Histogram {
   static const int nPer = 2 + 'H'-'B'; 

   TString name, title; 
   int nDim = 0;
   vector<TH1 *> hist;

};

struct HistoManager {



};

#define SF TString::Format 





void matching::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}


void matching::Histos::Init(TList *fOutput_)
{
    fOutput = fOutput_;

    TH1::SetDefaultSumw2();

    const int nAsym = 120;
    vector<double> asymBins;
    for(int i = 0; i <= nAsym; ++i) asymBins.push_back( 0.6 + 0.8 * i /(nAsym+0.0));


    vector<double> jecBins;
    for(double s = 0.9; s <= 1.6; s += 0.7/70)
        jecBins.push_back(s);

    vector<double> rhoBins;
    for(double s = 10.5; s <= 70; ++s)
        rhoBins.push_back(s);

    //Eta binning
    vector<double> etaBins;
    for(auto it = etaBinsPos.rbegin(); it != etaBinsPos.rend();  ++it)
        etaBins.push_back(-*it);
    etaBins.insert(etaBins.end(), etaBinsPos.begin()+1, etaBinsPos.end());
        
    for(auto x : etaBins)
        cout << x << endl;
    //exit(0);

    for(int i = 0; i < nPer; ++i) {
        char per = 'A' + i;
        for(int j = 0; j < 5; ++j) {//j=0 everything, others binned
            hBalEtaPt[j][i] = new TH3D(SF("hBalEtaPt%d_%c",j, per), SF("hBalEtaPt%d_%c",j, per), etaBins.size()-1, etaBins.data(), Ptbinning.size()-1, Ptbinning.data(),  asymBins.size()-1, asymBins.data()); //,  60, 0.8, 1.2);
            fOutput->Add(hBalEtaPt[j][i]);
        }
        hJetPt[i] = new TH1D(SF("hJetPt_%c",per), SF("hJetPt_%c",per), Ptbinning.size()-1, Ptbinning.data());
        fOutput->Add(hJetPt[i]);

        hJetPtInc[i] = new TH1D(SF("hJetPtInc_%c",per), SF("hJetPtInc_%c",per), Ptbinning.size()-1, Ptbinning.data());
        fOutput->Add(hJetPtInc[i]);


        for(int k = 0; k < 14; ++k) {
            hProf[k][i] = new TProfile2D(SF("hProf%d_%c",k+1, per), SF("hProf%d_%c",k+1, per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data());
            fOutput->Add(hProf[k][i]);
            hProfBB[k][i] = new TProfile2D(SF("hProfBB%d_%c",k+1, per), SF("hProfBB%d_%c",k+1, per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data());
            fOutput->Add(hProfBB[k][i]);
        }


        hRhopuppi[i] = new TH3D(SF("hRhopuppi_%c",per), SF("hRhopuppi_%c",per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data(),
                                                                          rhoBins.size()-1, rhoBins.data());
        fOutput->Add(hRhopuppi[i]);


        hRhochs[i] = new TH3D(SF("hRhochs_%c",per), SF("hRhochs_%c",per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data(),
                                                                          rhoBins.size()-1, rhoBins.data());
        fOutput->Add(hRhochs[i]);


        hEtaPtCHS[i] =  new TH2D(SF("hEtaPtCHS_%c",per), SF("hEtaPtCHS_%c",per), etaBins.size()-1, etaBins.data(),Ptbinning.size()-1, Ptbinning.data() );
        fOutput->Add(hEtaPtCHS[i]);
        hEtaPtPUPPI[i] = new TH2D(SF("hEtaPtPUPPI_%c",per), SF("hEtaPtPUPPI_%c",per), etaBins.size()-1, etaBins.data(),Ptbinning.size()-1, Ptbinning.data() );
        fOutput->Add(hEtaPtPUPPI[i]);
        hEtaPtPUPPIalone[i] = new TH2D(SF("hEtaPtPUPPIalone_%c",per), SF("hEtaPtPUPPIalone_%c",per), etaBins.size()-1, etaBins.data(),Ptbinning.size()-1, Ptbinning.data() );
        fOutput->Add(hEtaPtPUPPIalone[i]);

    }

}

void matching::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   h.Init(fOutput);
   lum.LoadLumis();


}



void CorrectJets(TTreeReaderValue<std::vector<QCDjet> > &Jets, double rho, JECs &jetCorrs ) {
    for (unsigned i = 0; i < Jets->size(); ++i) {
        double pt  = Jets->at(i).p4.Pt();
        double eta = Jets->at(i).p4.Eta();
        double area = Jets->at(i).area;

        //cout <<"R "<<  i <<" "<< pt << " "<< eta <<" "<< area <<" "<< rho <<  endl;
        double CorFactorRes, Unc;
        vector<string> dumy;

        double corr = jetCorrs.JEC_CHScorrections(pt, eta, area, rho, dumy, CorFactorRes, Unc);
        //cout << "Corr " << corr << endl;

        corr *= isLeptonVetoJet(Jets->at(i));
        if(withoutRes) corr /= CorFactorRes;

        Jets->at(i).p4.SetPt(pt * corr);

    }
    sort(Jets->begin(), Jets->end(), [](const QCDjet &a, const QCDjet &b) {return a.p4.Pt() > b.p4.Pt();});
    auto pEnd = std::remove_if (Jets->begin(), Jets->end(), [](const QCDjet &a) {return  a.p4.Pt() < 20;}); 
    Jets->resize(pEnd - Jets->begin());
}


void matching::DoMikkoMatching(int fileId)
{

    if(chsJets->size() < 2 || puppiJets->size() < 1) return;

    auto   CHStag   = chsJets->at(0).p4;
    auto   CHSprobe = chsJets->at(1).p4;


    double dPhi = abs(CHStag.Phi() - CHSprobe.Phi());
    dPhi = min(dPhi, 2*M_PI - dPhi);
    if (dPhi < 2.8) return;


    if(gRandom->Uniform() < 0.5) swap(CHStag, CHSprobe);


    if(abs(CHStag.Eta()) > 1.3) return;

    //Check trigger efficiency

    double wgt, wgtTot;
    int id;
    tie(wgt,wgtTot,id) = lum.GetWeightID(fileId, CHStag.Pt() );
    if(wgt == 0 || id < 0)  return;

    if(triggerBit->at(id) != 1) return;

    //Calculate the variable ala 
    //https://indico.cern.ch/event/544604/contributions/2210109/attachments/1294031/1928613/Giugno-17-2016_-_CaloScouting.pdf


    int m = -1;
    for(int j = 0; j < (int)puppiJets->size(); ++j) {
        double d2 = dist2(puppiJets->at(j).p4, CHSprobe);
        if(d2 < 0.2*0.2) {
            m = j;
            break;
        }
    }
    if(m == -1) return; //not found
    auto PUPPIprobe = puppiJets->at(m).p4;

    m = -1;
    for(int j = 0; j < (int)puppiJets->size(); ++j) {
        double d2 = dist2(puppiJets->at(j).p4, CHStag);
        if(d2 < 0.2*0.2) {
            m = j;
            break;
        }
    }
    auto PUPPItag = PUPPIprobe;
    if(m != -1) 
        PUPPItag = puppiJets->at(m).p4;




    double var1 = PUPPIprobe.Pt() /  CHStag.Pt();
    double var2 = CHSprobe.Pt() /  CHStag.Pt();
    double var3 = PUPPIprobe.Pt() / CHSprobe.Pt();
    double var4 = PUPPIprobe.Pt();

    double avgProbe = sqrt(CHSprobe.Pt() * PUPPIprobe.Pt());

    double pTtag = CHStag.Pt();
    double eta   = PUPPIprobe.Eta();
    auto fillHist = [&](array<TProfile2D*,nPer> *hProf,  int fId, double w) {
        //cout << "Fill begin " << var1 << " "<< var2 <<" "<< var3 <<" "<<var4 << endl;
        //cout << "b"<<endl;
        hProf[0][fId]->Fill(eta, pTtag, var1, w);
        hProf[1][fId]->Fill(eta, pTtag, var2, w);
        hProf[2][fId]->Fill(eta, pTtag, var3, w);
        hProf[3][fId]->Fill(eta, pTtag, var4, w);
        hProf[4][fId]->Fill(eta, PUPPIprobe.Pt(), var3, w); //resambles the old way
        hProf[5][fId]->Fill(eta, CHSprobe.Pt(),   var3, w); 
        hProf[6][fId]->Fill(eta, avgProbe,        var3, w); 
        hProf[7][fId]->Fill(CHSprobe.Eta(), pTtag,   CHSprobe.Pt(), w); 

        hProf[8][fId]->Fill(eta, PUPPIprobe.Pt(), CHSprobe.Pt(), w); //resambles the old way
        hProf[9][fId]->Fill(eta, CHSprobe.Pt(),   PUPPIprobe.Pt(), w); 

        if(m != -1) { //puppiTag exist
            double pTtagP = PUPPItag.Pt(); 
            hProf[10][fId]->Fill(eta, pTtagP, var3, w);
            hProf[11][fId]->Fill(eta, pTtagP, var4, w);
            hProf[12][fId]->Fill(CHSprobe.Eta(), pTtagP,   CHSprobe.Pt(), w); 
        }

    };

    fillHist(h.hProf, 0, wgtTot);
    fillHist(h.hProf, fileId, wgt);

    bool isBB = true;
    if(chsJets->size() >= 3) {
        double pM = (chsJets->at(0).p4.Pt() + chsJets->at(1).p4.Pt())/2;
        double p3 =  chsJets->at(2).p4.Pt();
        if(p3 / pM > 0.3) //skip alpha > 0.3
            isBB = false;
    }

    if(isBB) {
        fillHist(h.hProfBB, 0, wgtTot);
        fillHist(h.hProfBB, fileId, wgt);
    }


}



Bool_t matching::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   //cout << "I am reading event " << entry << endl;

   fReader.SetEntry(entry);

    //cout << "Init start " <<entry <<" "<< __LINE__ << endl;
   if(entry % 100000 == 0)
       cout << "Event " << entry << endl;



    TString fileName =  fReader.GetTree()->GetCurrentFile()->GetName();
    int fileId = fileName[fileName.Length()-6] -'A';
    if(fileName != currFile) {
        cout << "Filename is " <<   fReader.GetTree()->GetCurrentFile()->GetName() << endl;

        //string jecTag = "Summer16_07Aug2017";
        //int version = 4;
        char period =  fileName[fileName.Length()-6];
        vector<string> dumy;
        bool isMC = false;

        jetCorrsCHS.Init(isMC, jecTagCHS, period, versionCHS, "AK4PFchs", "", dumy);
        jetCorrsPUPPI.Init(isMC, jecTagPUPPI, period, versionPUPPI, "AK4PFPuppi", "", dumy);
        currFile = fileName;
        cout << "JEC init succesfull" << endl;
    }

    //cout << "Init start " <<entry <<" "<< __LINE__ << endl;

    CorrectJets(chsJets,   *rho, jetCorrsCHS);
    CorrectJets(puppiJets, *rho, jetCorrsPUPPI);


    const double jetSel = 5;


    if(chsJets->size() < 1 || puppiJets->size() < 1) return 0;


    DoMikkoMatching(fileId);




    //cout << "Init start " <<entry <<" "<< __LINE__ << endl;

    double wgt, wgtTot;
    int id;
    tie(wgt,wgtTot,id) = lum.GetWeightID(fileId, chsJets->at(0).p4.Pt() );
    if(wgt == 0 || id < 0)  return 0;
    h.w = wgt;
    h.wTot = wgtTot;
    h.fileId = fileId;

    //cout << "Info " << CHSjetPt[0]<<" "<< id <<" "<< (*CHStriggerBit).at(id) << endl;
    if(triggerBit->at(id) != 1) return 0;

    h.Fill1D(h.hJetPt, chsJets->at(0).p4.Pt());


    for(unsigned i = 0; i < chsJets->size(); ++i) {
        h.Fill1D(h.hJetPtInc, chsJets->at(i).p4.Pt());
    }


    //cout << "Init start " <<entry <<" "<< __LINE__ << endl;
    //cout << "Second part OK" << endl;



    //cout << CHSjetPt[0] << " "<< CHSjetPt[0]*(1+CHSjetUnc[0]) <<" "<< CHSjetPt[0]*(1-CHSjetUnc[0]) << endl;


    for (unsigned i = 0; i < puppiJets->size(); ++i) {
        double eta = puppiJets->at(i).p4.Eta();
        double pt = puppiJets->at(i).p4.Pt();
        double area = puppiJets->at(i).area;

        //double CorFactorRes, Unc;
        //vector<string> dumy;
        //double corr = jetCorrsPUPPI.JEC_CHScorrections(pt, eta, area, *rho, dumy, CorFactorRes, Unc);

        h.Fill3D(h.hRhopuppi, eta, pt, *rho);
    }
    //cout << "Init start " <<entry <<" "<< __LINE__ << endl;
    for (unsigned i = 0; i < chsJets->size(); ++i) {
        double eta = chsJets->at(i).p4.Eta();
        double pt = chsJets->at(i).p4.Pt();
        double area = chsJets->at(i).area;

        //double CorFactorRes, Unc;
        //vector<string> dumy;
        //double corr = jetCorrsCHS.JEC_CHScorrections(pt, eta, area, *rho, dumy, CorFactorRes, Unc);

        h.Fill3D(h.hRhochs, eta, pt, *rho);

        if (pt >= jetSel) {
            h.Fill2D(h.hEtaPtCHS, eta, pt);
        }
    }
    //cout << "Init start " <<entry <<" "<< __LINE__ << endl;

    set<int> Indx;
    for(unsigned i = 0; i < chsJets->size(); ++i)
        Indx.insert(i);


    for(unsigned i = 0; i < puppiJets->size(); ++i) {
        double pt   = puppiJets->at(i).p4.Pt();
        double Eta = puppiJets->at(i).p4.Eta();
        if(pt < jetSel) continue;

        h.Fill2D(h.hEtaPtPUPPI, Eta, pt);

        /*
           int m = -1;
           for(int j = 0; j < CHSjetPt.GetSize(); ++j) {
           double d2 = dist2(PUPPIjetEta[i], PUPPIjetPhi[i], CHSjetEta[j], CHSjetPhi[j]);
           if(d2 < 0.2*0.2) {
           m = j;
           break;
           }
           }
           */

        int m = -1;
        for(int ind :  Indx) {
            //cout << ind << endl;
            double d2 = dist2(puppiJets->at(i).p4, chsJets->at(ind).p4);
            if(d2 < 0.2*0.2) {
                m = ind;
                Indx.erase(m);
                break;
            }
        }


        if(m != -1) {
            //cout << "match " << PUPPIjetPt[i] << " "<< PUPPIjetPt[i] / CHSjetPt[m] << endl;
            double r = puppiJets->at(i).p4.Pt() / chsJets->at(m).p4.Pt();
            h.Fill3D(h.hBalEtaPt[0], Eta, pt, r);
            if(*rho < 15)
                h.Fill3D(h.hBalEtaPt[1], Eta, pt, r);
            else if(*rho < 22)
                h.Fill3D(h.hBalEtaPt[2], Eta, pt, r);
            else if(*rho < 30)
                h.Fill3D(h.hBalEtaPt[3], Eta, pt, r);
            else
                h.Fill3D(h.hBalEtaPt[4], Eta, pt, r);

            //cout <<"Event " <<  PUPPIjetJEC[i] << " "<< CHSjetJEC[m] << endl;
        }
        else { //not matched
            h.Fill2D(h.hEtaPtPUPPIalone, Eta, pt);
        }
    }


   return kTRUE;
}

void matching::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void matching::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.


    cout << "Saving histos" << endl;
    //return;

    gRandom->SetSeed(0);
    TString outputFilename = SF("histos/%sV%d__%sV%d%s_%d.root", jecTagCHS.c_str(), versionCHS,    jecTagPUPPI.c_str(), versionPUPPI,
           withoutRes ? "noRes" : "", gRandom->Integer(999999));
    TFile* outputFile = new TFile(outputFilename,"recreate");

    for(const auto&& obj: *fOutput) {
        if (obj->InheritsFrom("TH1"))
            obj->Write();
    }

    outputFile->Close();

}
