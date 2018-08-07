#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TRandom.h>
#include <vector>

#include "utils.h"
#include "JEC/treeProducer/interface/QCDjet.h"
#include "JEC.h"
#include "lumi.h"


using namespace std::experimental::filesystem;
using namespace std;





JECs * GetDataJEC (int runNo)
{
    if (runNo >= 304911) {
        static JECs * jecs2017F = new JECs("Fall17_17Nov2017F_V6_DATA");
        return jecs2017F;
    }
    else if (runNo >= 303435) { 
        static JECs * jecs2017E = new JECs("Fall17_17Nov2017E_V6_DATA");
        return jecs2017E;
    }
    else if (runNo >= 302030) { 
        static JECs * jecs2017D = new JECs("Fall17_17Nov2017D_V6_DATA");
        return jecs2017D;
    }
    else if (runNo >= 299337) { 
        static JECs * jecs2017C = new JECs("Fall17_17Nov2017C_V6_DATA");
        return jecs2017C;
    }
    else if (runNo >= /*297179*/297000) { // TODO: CHECK!
        static JECs * jecs2017B = new JECs("Fall17_17Nov2017B_V6_DATA");
        return jecs2017B;
    }
    else if (runNo >= 294034) {
        cout << runNo << " belongs to Run2017A, which is not a physics run!" << endl;
        exit(EXIT_FAILURE);
    }
    else if (runNo >= 278802) {
        static JECs * jecs2016GH = new JECs("Summer16_07Aug2017GH_V11_DATA");
        return jecs2016GH;
    }
    else if (runNo >= 276831) {
        static JECs * jecs2016EF = new JECs("Summer16_07Aug2017EF_V11_DATA");
        return jecs2016EF;
    }
    else if (runNo >= 272007) {
        static JECs * jecs2016BCD = new JECs("Summer16_07Aug2017BCD_V11_DATA");
        return jecs2016BCD;
    }
    else 
        exit(EXIT_FAILURE);
}

void ApplyJEC (int runNo, double rho, vector<QCDjet> & jets)
{
    JECs * jecs = GetDataJEC(runNo); 
    //cout << runNo << '\t' << oldera << '\t' << era << '\t' << endl;
    //cout << runNo << '\t' << rho << '\t' << pt << '\t' << eta << '\t' << corr << '\t';

    for(auto &jet : jets) {
        double corr = jecs->getJEC(jet, rho);
        jet.p4.Scale(corr);
    }
    sort(jets.begin(), jets.end(), [](QCDjet &a, QCDjet &b) {return a.p4.Pt() > b.p4.Pt();});
}



void procesor (TString input, TString output, int nSplit = 1, int nNow = 0)
{
    //gROOT->SetBatch();
    //Get old file, old tree and set top branch address
    TChain * inputChain = new TChain("ak4/events"); // arg = path inside the ntuple
    for(auto f : GetFiles(input))
        if(!f.Contains("failed")) inputChain->Add(f);

    int runNo;
    float rho;
    vector<bool> *triggerBit=nullptr;
    vector<QCDjet> *chsJets=nullptr, *puppiJets=nullptr;
    inputChain->SetBranchAddress("runNo", &runNo);
    inputChain->SetBranchAddress("rho", &rho);
    inputChain->SetBranchAddress("triggerBit", &triggerBit);
    inputChain->SetBranchAddress("chsJets", &chsJets);
    inputChain->SetBranchAddress("puppiJets", &puppiJets);


    // saving MC weight 
    TFile *outputFile = TFile::Open(output, "RECREATE");
    TTree *outputTree = inputChain->CloneTree(0);
    float wgt_, wgtTot_;
    outputTree->Branch("wgt"    ,&wgt_      ,"wgt/F");
    outputTree->Branch("wgtTot"    ,&wgtTot_   ,"wgtTot/F");

    cout << "Calculating #entries" << endl;
    Long64_t N = inputChain->GetEntries();
    cout << input << '\t' << N << '\t' << output << endl;

    Luminosity lum;
    lum.LoadLumis();

    auto range = splitRange(N, nSplit, nNow);

    for (int i=range.first; i < range.second; ++i) { //Event LOOP --start
        PrintCounterAndIncrement(output, range);
        inputChain->GetEntry(i);

        //JEC corrections
        ApplyJEC(runNo, rho, *chsJets);
        ApplyJEC(runNo, rho, *puppiJets);

        //For use in probe+tag method at least two jets are needed
        if(chsJets->size() < 2) continue;
        //Randomly swap leading and subleading
        if(gRandom->Uniform() < 0.5) swap(chsJets->at(0), chsJets->at(1));
        //cout <<"Radek " <<  chsJets->at(0).p4.Pt() <<" "<<  chsJets->at(1).p4.Pt() << endl;

        int perID = getPer(runNo) - 'A';
        int id;
        double wgt, wgtTot;
        tie(wgt,wgtTot,id) = lum.GetWeightID(perID, chsJets->at(0).p4.Pt()); //check that tag jet is within trigger acceptance
        //cout << wgt << " "<< wgtTot << endl;
        if(wgt == 0 || id < 0)  continue;
        if(triggerBit->at(id) != 1) continue;
        wgt_ = wgt; wgtTot_ = wgtTot;

        outputTree->Fill();
    } //Event LOOP --start

    /* closing */
    outputTree->AutoSave();
    outputFile->Close();
    delete inputChain;
    cout << "Done " << nNow << endl;
}

int main (int argc, char * argv[])
{
    if      (argc == 3) procesor(argv[1], argv[2]);
    else if (argc == 4) procesor(argv[1], argv[2], atoi(argv[3]));
    else if (argc == 5) procesor(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
    else {
        cout << "Wrong number of arguments" << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
