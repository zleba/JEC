#ifndef _filler_h
#define _filler_h

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile2D.h"


static const vector<double> etaBinsPos =  {0., 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172,  2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.538, 4.716, 4.889, 5.191};


static const vector<double> Ptbinning = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638,     686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,  2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961,     5220, 5492, 5777, 6076, 6389, 6717,    7000};


    static const int nPer = 8;
    struct Histos {
        //TList  *fOutput;
        double w, wTot;
        int fileId;
        array<TH3D*,nPer> hBalEtaPt[5];
        array<TH1D*,nPer> hJetPt;
        array<TH1D*,nPer> hJetPtInc;

        array<TProfile2D*,nPer> hProf[14], hProfBB[14];
                                          
        array<TH3D*,nPer> hRhopuppi;
        array<TH3D*,nPer> hRhochs;
                                          
        array<TH2D*,nPer> hEtaPtCHS;
        array<TH2D*,nPer> hEtaPtPUPPI;
        array<TH2D*,nPer> hEtaPtPUPPIalone;
        void Init();

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


void Histos::Init()
{

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
            hBalEtaPt[j][i] = new TH3D(Form("hBalEtaPt%d_%c",j, per), Form("hBalEtaPt%d_%c",j, per), etaBins.size()-1, etaBins.data(), Ptbinning.size()-1, Ptbinning.data(),  asymBins.size()-1, asymBins.data()); //,  60, 0.8, 1.2);
        }
        hJetPt[i] = new TH1D(Form("hJetPt_%c",per), Form("hJetPt_%c",per), Ptbinning.size()-1, Ptbinning.data());

        hJetPtInc[i] = new TH1D(Form("hJetPtInc_%c",per), Form("hJetPtInc_%c",per), Ptbinning.size()-1, Ptbinning.data());


        for(int k = 0; k < 14; ++k) {
            hProf[k][i] = new TProfile2D(Form("hProf%d_%c",k+1, per), Form("hProf%d_%c",k+1, per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data());
            hProfBB[k][i] = new TProfile2D(Form("hProfBB%d_%c",k+1, per), Form("hProfBB%d_%c",k+1, per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data());
        }

        hRhopuppi[i] = new TH3D(Form("hRhopuppi_%c",per), Form("hRhopuppi_%c",per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data(),
                                                                          rhoBins.size()-1, rhoBins.data());


        hRhochs[i] = new TH3D(Form("hRhochs_%c",per), Form("hRhochs_%c",per), etaBins.size()-1, etaBins.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data(),
                                                                          rhoBins.size()-1, rhoBins.data());

        hEtaPtCHS[i] =  new TH2D(Form("hEtaPtCHS_%c",per), Form("hEtaPtCHS_%c",per), etaBins.size()-1, etaBins.data(),Ptbinning.size()-1, Ptbinning.data() );
        hEtaPtPUPPI[i] = new TH2D(Form("hEtaPtPUPPI_%c",per), Form("hEtaPtPUPPI_%c",per), etaBins.size()-1, etaBins.data(),Ptbinning.size()-1, Ptbinning.data() );
        hEtaPtPUPPIalone[i] = new TH2D(Form("hEtaPtPUPPIalone_%c",per), Form("hEtaPtPUPPIalone_%c",per), etaBins.size()-1, etaBins.data(),Ptbinning.size()-1, Ptbinning.data() );

    }

}

#endif
