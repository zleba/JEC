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





JECs * GetDataJEC (int runNo, string jetType)
{
    char per = getPer(runNo);

    static map<char,JECs*> jecMapCHS4, jecMapCHS8, jecMapPUPPI4;

    if (jecMapCHS4.count(per) == 0) {
        if ('B' <= per  && per <= 'D') {
            jecMapCHS4[per]   = new JECs("Summer16_07Aug2017BCD_V11_DATA", "AK4PFchs");
            jecMapCHS8[per]   = new JECs("Summer16_07Aug2017BCD_V11_DATA", "AK8PFchs");
            jecMapPUPPI4[per] = new JECs("Summer16_07Aug2017BCD_V11_DATA", "AK4PFPuppi");
        }
        else if ('E' <= per  && per <= 'F') {
            jecMapCHS4[per]   = new JECs("Summer16_07Aug2017EF_V11_DATA", "AK4PFchs");
            jecMapCHS8[per]   = new JECs("Summer16_07Aug2017EF_V11_DATA", "AK8PFchs");
            jecMapPUPPI4[per] = new JECs("Summer16_07Aug2017EF_V11_DATA", "AK4PFPuppi");
        }
        else if('G' <= per  && per <= 'H') {
            jecMapCHS4[per]   = new JECs("Summer16_07Aug2017GH_V11_DATA", "AK4PFchs");
            jecMapCHS8[per]   = new JECs("Summer16_07Aug2017GH_V11_DATA", "AK8PFchs");
            jecMapPUPPI4[per] = new JECs("Summer16_07Aug2017GH_V11_DATA", "AK4PFPuppi");
        }
        else {
            cout << "Wrong period "<< per  << endl;
            exit(1);
        }
    }

    if (jetType == "AK4PFchs")
        return jecMapCHS4.at(per);
    else if (jetType == "AK8PFchs")
        return jecMapCHS8.at(per);
    else if (jetType == "AK4PFPuppi")
        return jecMapPUPPI4.at(per);
    else {
        cout << "Wrong jet type" << endl;
        exit(1);
    }
}

void ApplyJEC (int runNo, double rho, vector<QCDjet> & jets, string jetType)
{
    JECs * jecs = GetDataJEC(runNo, jetType); 
    //cout << runNo << '\t' << oldera << '\t' << era << '\t' << endl;
    //cout << runNo << '\t' << rho << '\t' << pt << '\t' << eta << '\t' << corr << '\t';

    for(auto &jet : jets) {
        double corr = jecs->getJEC(jet, rho);
        //cout << "RADEK " << jetType << " "<< jet.p4.Pt() << " "<< corr << endl;
        jet.p4.Scale(corr);
    }
    sort(jets.begin(), jets.end(), [](QCDjet &a, QCDjet &b) {return a.p4.Pt() > b.p4.Pt();});
}



void procesor (TString input, TString output, int nSplit = 1, int nNow = 0)
{
    //gROOT->SetBatch();
    //Get old file, old tree and set top branch address
    TChain * inputChain = new TChain("alljets/events"); // arg = path inside the ntuple
    for(auto f : GetFiles(input))
        if(!f.Contains("failed")) inputChain->Add(f);

    int runNo;
    float rho;
    vector<bool> *triggerBit=nullptr;
    vector<QCDjet> *chs4Jets=nullptr, *chs8Jets=nullptr, *puppi4Jets=nullptr;
    inputChain->SetBranchAddress("runNo", &runNo);
    inputChain->SetBranchAddress("rho", &rho);
    inputChain->SetBranchAddress("triggerBit", &triggerBit);
    inputChain->SetBranchAddress("chs4Jets", &chs4Jets);
    inputChain->SetBranchAddress("chs8Jets", &chs8Jets);
    inputChain->SetBranchAddress("puppi4Jets", &puppi4Jets);


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

        //if(chs8Jets->size() > 0) cout << "TEST1 " << chs8Jets->at(0).p4.Pt() << endl;

        //JEC corrections
        ApplyJEC(runNo, rho, *chs4Jets,   "AK4PFchs");
        ApplyJEC(runNo, rho, *chs8Jets,   "AK8PFchs");
        ApplyJEC(runNo, rho, *puppi4Jets, "AK4PFPuppi");

        //if(chs8Jets->size() > 0) cout << "TEST2 " << chs8Jets->at(0).p4.Pt() << endl;

        //For use in probe+tag method at least two jets are needed
        if(chs4Jets->size() < 2) continue;
        //Randomly swap leading and subleading
        if(gRandom->Uniform() < 0.5) swap(chs4Jets->at(0), chs4Jets->at(1));
        //cout <<"Radek " <<  chs4Jets->at(0).p4.Pt() <<" "<<  chs4Jets->at(1).p4.Pt() << endl;

        int perID = getPer(runNo) - 'A';
        int id;
        double wgt, wgtTot;
        tie(wgt,wgtTot,id) = lum.GetWeightID(perID, chs4Jets->at(0).p4.Pt()); //check that tag jet is within trigger acceptance
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
    cout << "Done " << nSplit <<" "<< nNow << endl;
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
