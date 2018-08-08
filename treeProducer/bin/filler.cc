#include "utils.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TRandom.h"
#include <set>
#include "JEC/treeProducer/interface/QCDjet.h"
#include "filler.h"

using namespace std;



//double dist2(double eta1, double phi1, double eta2, double phi2)
double dist2(const ROOT::Math::PtEtaPhiM4D<float> &j1, const ROOT::Math::PtEtaPhiM4D<float> &j2)
{
    double dPhi = abs(j2.Phi() - j1.Phi());
    dPhi = min(dPhi, 2*M_PI - dPhi);
    double dEta = abs(j2.Eta() - j1.Eta());

    return (dEta*dEta + dPhi*dPhi);
}

Histos h;


void DoMikkoMatching(int runNo, vector<QCDjet> *chsJets, vector<QCDjet> *testJets, double wgt, double wgtTot)
{
    int fileId =  getPer(runNo) - 'A';

    if(chsJets->size() < 2 || testJets->size() < 1) return;

    auto   CHStag   = chsJets->at(0).p4;
    auto   CHSprobe = chsJets->at(1).p4;


    double dPhi = abs(CHStag.Phi() - CHSprobe.Phi());
    dPhi = min(dPhi, 2*M_PI - dPhi);
    if (dPhi < 2.8) return;


    if(abs(CHStag.Eta()) > 1.3) return;


    //Calculate the variable ala 
    //https://indico.cern.ch/event/544604/contributions/2210109/attachments/1294031/1928613/Giugno-17-2016_-_CaloScouting.pdf

    int m = -1;
    for(int j = 0; j < (int) testJets->size(); ++j) {
        double d2 = dist2(testJets->at(j).p4, CHSprobe);
        if(d2 < 0.2*0.2) {
            m = j;
            break;
        }
    }
    if(m == -1) return; //not found
    auto TESTprobe = testJets->at(m).p4;

    m = -1;
    for(int j = 0; j < (int)testJets->size(); ++j) {
        double d2 = dist2(testJets->at(j).p4, CHStag);
        if(d2 < 0.2*0.2) {
            m = j;
            break;
        }
    }
    auto TESTtag = TESTprobe;
    if(m != -1) TESTtag = testJets->at(m).p4;



    double var1 = TESTprobe.Pt() /  CHStag.Pt();
    double var2 = CHSprobe.Pt() /  CHStag.Pt();
    double var3 = TESTprobe.Pt() / CHSprobe.Pt();
    double var4 = TESTprobe.Pt();

    double avgProbe = sqrt(CHSprobe.Pt() * TESTprobe.Pt());

    double pTtag = CHStag.Pt();
    double eta   = TESTprobe.Eta();
    auto fillHist = [&](array<TProfile2D*,nPer> *hProf,  int fId, double w) {
        //cout << "Fill begin " << var1 << " "<< var2 <<" "<< var3 <<" "<<var4 << endl;
        //cout << "b"<<endl;
        hProf[0][fId]->Fill(eta, pTtag, var1, w);
        hProf[1][fId]->Fill(eta, pTtag, var2, w);
        hProf[2][fId]->Fill(eta, pTtag, var3, w);
        hProf[3][fId]->Fill(eta, pTtag, var4, w);
        hProf[4][fId]->Fill(eta, TESTprobe.Pt(), var3, w); //resambles the old way
        hProf[5][fId]->Fill(eta, CHSprobe.Pt(),   var3, w); 
        hProf[6][fId]->Fill(eta, avgProbe,        var3, w); 
        hProf[7][fId]->Fill(CHSprobe.Eta(), pTtag,   CHSprobe.Pt(), w); 

        hProf[8][fId]->Fill(eta, TESTprobe.Pt(), CHSprobe.Pt(), w); //resambles the old way
        hProf[9][fId]->Fill(eta, CHSprobe.Pt(),   TESTprobe.Pt(), w); 

        if(m != -1) { //puppiTag exist
            double pTtagP = TESTtag.Pt(); 
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




void filler(TString input, TString output, int nSplit = 1, int nNow = 0)
{


    TChain * oldchain = new TChain("events"); // arg = path inside the ntuple
    oldchain->Add(input);


    int runNo;
    float rho;
    float wgt, wgtTot;
    vector<QCDjet> *chsJets=nullptr, *testJets=nullptr;
    oldchain->SetBranchAddress("runNo", &runNo);
    oldchain->SetBranchAddress("rho", &rho);
    oldchain->SetBranchAddress("wgt", &wgt);
    oldchain->SetBranchAddress("wgtTot", &wgtTot);
    //oldchain->SetBranchAddress("triggerBit", &triggerBit);
    oldchain->SetBranchAddress("chs4Jets", &chsJets);
    oldchain->SetBranchAddress("chs8Jets", &testJets);

    TFile *newfile = TFile::Open(output, "RECREATE");
    h.Init();

    cout << h.hProf[0][0] << endl;

   Long64_t N = oldchain->GetEntries();
   auto range = splitRange(N, nSplit, nNow);
   for (int i=range.first; i < range.second; ++i) {
        PrintCounterAndIncrement(output, range);
        oldchain->GetEntry(i);

        const double jetSel = 5;

        if(chsJets->size() < 1 || testJets->size() < 1) continue;


        DoMikkoMatching(runNo, chsJets, testJets, wgt, wgtTot);

        h.w = wgt;
        h.wTot = wgtTot;

        h.Fill1D(h.hJetPt, chsJets->at(0).p4.Pt());


        for(unsigned i = 0; i < chsJets->size(); ++i) {
            h.Fill1D(h.hJetPtInc, chsJets->at(i).p4.Pt());
        }


        for (unsigned i = 0; i < testJets->size(); ++i) {
            double eta = testJets->at(i).p4.Eta();
            double pt = testJets->at(i).p4.Pt();

            //double CorFactorRes, Unc;
            //vector<string> dumy;
            //double corr = jetCorrsPUPPI.JEC_CHScorrections(pt, eta, area, *rho, dumy, CorFactorRes, Unc);

            h.Fill3D(h.hRhopuppi, eta, pt, rho);
        }
        //cout << "Init start " <<entry <<" "<< __LINE__ << endl;


        set<int> Indx;
        for(unsigned i = 0; i < testJets->size(); ++i)
            Indx.insert(i);

        for (unsigned i = 0; i < chsJets->size(); ++i) {
            double eta = chsJets->at(i).p4.Eta();
            double pt = chsJets->at(i).p4.Pt();
            //double area = chsJets->at(i).area;

            //double CorFactorRes, Unc;
            //vector<string> dumy;
            //double corr = jetCorrsCHS.JEC_CHScorrections(pt, eta, area, *rho, dumy, CorFactorRes, Unc);

            h.Fill3D(h.hRhochs, eta, pt, rho);

            if (pt >= jetSel) {
                h.Fill2D(h.hEtaPtCHS, eta, pt);

                int m = -1;
                for(int ind :  Indx) {
                    //cout << ind << endl;
                    double d2 = dist2(chsJets->at(i).p4, testJets->at(ind).p4);
                    if(d2 < 0.2*0.2) {
                        m = ind;
                        Indx.erase(m);
                        break;
                    }
                }
                //cout << "CHS jet " << pt << " "<< m << endl;
                if(m == -1) {
                    h.Fill2D(h.hEtaPtCHSalone, eta, pt);
                }
            }
        }
        //cout << "Init start " <<entry <<" "<< __LINE__ << endl;

        Indx.clear();
        for(unsigned i = 0; i < chsJets->size(); ++i)
            Indx.insert(i);


        //Loop over tested jets and check mathing properties
        for(unsigned i = 0; i < testJets->size(); ++i) {
            double pt   = testJets->at(i).p4.Pt();
            double eta = testJets->at(i).p4.Eta();
            if(pt < jetSel) continue;

            //cout << "Filling pupi" << 
            h.Fill2D(h.hEtaPtPUPPI, eta, pt);

            int m = -1;
            for(int ind :  Indx) {
                //cout << ind << endl;
                double d2 = dist2(testJets->at(i).p4, chsJets->at(ind).p4);
                if(d2 < 0.2*0.2) {
                    m = ind;
                    Indx.erase(m);
                    break;
                }
            }

            if(m != -1) {
                //cout << "match " << PUPPIjetPt[i] << " "<< PUPPIjetPt[i] / CHSjetPt[m] << endl;
                double r = testJets->at(i).p4.Pt() / chsJets->at(m).p4.Pt();
                h.Fill3D(h.hBalEtaPt[0], eta, pt, r);
                if(rho < 15)
                    h.Fill3D(h.hBalEtaPt[1], eta, pt, r);
                else if(rho < 22)
                    h.Fill3D(h.hBalEtaPt[2], eta, pt, r);
                else if(rho < 30)
                    h.Fill3D(h.hBalEtaPt[3], eta, pt, r);
                else
                    h.Fill3D(h.hBalEtaPt[4], eta, pt, r);

                //cout <<"Event " <<  PUPPIjetJEC[i] << " "<< CHSjetJEC[m] << endl;
            }
            else { //not matched
                h.Fill2D(h.hEtaPtPUPPIalone, eta, pt);
            }
        }

    }

    newfile->Write();
    newfile->Close();

}


int main (int argc, char * argv[])
{
    if      (argc == 3) filler(argv[1], argv[2]);
    else if (argc == 4) filler(argv[1], argv[2], atoi(argv[3]));
    else if (argc == 5) filler(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
    else {
        cout << "Wrong number of arguments" << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
