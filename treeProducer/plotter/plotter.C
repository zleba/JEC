#define SF TString::Format

#include "/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/libPlottingHelper.so)

//#include "utils/RemoveOverlaps.h"

#include "TH3D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TStyle.h"

#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include <algorithm>



using namespace PlottingHelper;


double polFun(double pt, double y, const double *x);


struct PLOTTER {

    const int nPer = 8;
    int perID;

    vector<TH3D *> hBalEtaPtAll[5];
    vector<TH3D *> hBalEtaPtAllUp;
    vector<TH3D *> hBalEtaPtAllDn;
    vector<TH1D *> hJetPtAll;

    vector<TH3D *> hJECpuppi;
    vector<TH3D *> hJECchs;

    vector<TH2D *> hEtaPtCHS;
    vector<TH2D *> hEtaPtCHSalone;
    vector<TH2D *> hEtaPtPUPPI;
    vector<TH2D *> hEtaPtPUPPIalone;

    vector<TProfile2D*> hProf[13], hProfBB[13];

    TF2 *f2;
    TH2D * hFit;

    TString outName = "newPlots.pdf";

    public:
    void Init(TString inFile, TString outFile) {
        outName = outFile;
        //TString fileName = "histoDir/jetsAll.root";
        TFile *fIn = TFile::Open(inFile, "READ");

        for(int j = 0; j < sizeof(hBalEtaPtAll)/sizeof(hBalEtaPtAll[0]); ++j)
            hBalEtaPtAll[j].resize(nPer);
        hBalEtaPtAllUp.resize(nPer);
        hBalEtaPtAllDn.resize(nPer);
        hJetPtAll.resize(nPer);
        hJECpuppi.resize(nPer);
        hJECchs.resize(nPer);

        hEtaPtCHS.resize(nPer);
        hEtaPtCHSalone.resize(nPer);
        hEtaPtPUPPI.resize(nPer);
        hEtaPtPUPPIalone.resize(nPer);

        for(auto & h : hProf) h.resize(nPer);
        for(auto & h : hProfBB) h.resize(nPer);

        for(int i = 0; i < nPer; ++i) {
            char per = 'A' + i;
            for(int j = 0; j < sizeof(hBalEtaPtAll)/sizeof(hBalEtaPtAll[0]); ++j) {
                hBalEtaPtAll[j][i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt%d_%c", j, per)));
                //hBalEtaPtAll[j][i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt_%c",  per)));
                assert(hBalEtaPtAll[j][i]);
            }

            for(int j = 0; j < sizeof(hProf)/sizeof(hProf[0]); ++j) {
                hProf[j][i] = dynamic_cast<TProfile2D*>(fIn->Get(SF("hProf%d_%c",j+1, per)));
                hProfBB[j][i] = dynamic_cast<TProfile2D*>(fIn->Get(SF("hProfBB%d_%c",j+1, per)));
                assert(hProf[j][i]);
                assert(hProfBB[j][i]);
            }


            hBalEtaPtAllUp[i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt0_%c",per)));
            assert(hBalEtaPtAllUp[i]);
            hBalEtaPtAllDn[i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt0_%c",per)));
            assert(hBalEtaPtAllDn[i]);

            hJetPtAll[i] = dynamic_cast<TH1D*>(fIn->Get(SF("hJetPt_%c",per)));
            assert(hJetPtAll[i]);

            hJECpuppi[i] = dynamic_cast<TH3D*>(fIn->Get(SF("hJECpuppi_%c",per)));
            assert(hJECpuppi[i]);
            hJECchs[i] = dynamic_cast<TH3D*>(fIn->Get(SF("hJECchs_%c",per)));
            assert(hJECchs[i]);


            hEtaPtCHS[i] = dynamic_cast<TH2D*>(fIn->Get(SF("hEtaPtCHS_%c",per)));
            assert(hEtaPtCHS[i]);
            hEtaPtCHSalone[i] = dynamic_cast<TH2D*>(fIn->Get(SF("hEtaPtCHSalone_%c",per)));
            assert(hEtaPtCHSalone[i]);

            hEtaPtPUPPI[i] = dynamic_cast<TH2D*>(fIn->Get(SF("hEtaPtPUPPI_%c",per)));
            assert(hEtaPtPUPPI[i]);
            hEtaPtPUPPIalone[i] = dynamic_cast<TH2D*>(fIn->Get(SF("hEtaPtPUPPIalone_%c",per)));
            assert(hEtaPtPUPPIalone[i]);


        }
    }

    void AsymmetryEtaDep();
    void AsymmetryPtDep();
    void AsymmetryEtaPtDep(int eta1, int eta2);
    void AsymmetryEtaPtTimeDep(int eta1, int eta2);
    void AsymmetryEtaPtPileUpDep(int eta1, int eta2);

    void PlotProfiles(int type, int per, bool doSigma);
    void PlotProfilesOverlap();
    pair<TH1D*,TH1D*> PlottMatchingCorr(int direction, bool doFit = true);

    void CompareSigmas();
    void PlotVar1overVar2();

    vector<double> MeanAsym(int shift=0, TString style="");
    void JEC();
    void PtEtaDep();
    void Unmatched(TString type = "tested");
    void MatchingFactorsTimeDep(int rho=0);
    pair<TH1D*,TH1D*>  MatchingFactorsPileUpDep(bool plotFit, TString plotType = "time");
    void  PlotFit(TH1D *h84, TH1D *h1000);
    void  PlotFitFour(TH1D *h84p, TH1D *h1000p, TH1D *h84n, TH1D *h1000n);
    void MatchingFactorsYbased(vector<double> res, bool allPer = true);

};


void PLOTTER::PlotProfiles(int type, int per, bool doSigma = false)
{
    if(type == 8 || type == 9)
        doSigma = true;


    const int pMax = per > 0 ? 8  : 1;
    TCanvas *can = new TCanvas("can", "can", 800, 500);

    if((type == 3 || type == 7) && !doSigma)
        SetLeftRight(0.28, 0.13);
    else
        SetLeftRight(0.08, 0.03);

    DividePad(vector<double>(7,1.), vector<double>(6,1.));

    gStyle->SetOptStat(0);

    for(int i = 1; i <= 7*6; ++i) {
        int ii = 7*6 + i;
        can->cd(i);
        TProfile *prof[8], *profBB[8];
        TProfile *Nprof[8], *NprofBB[8];
        for(int p = 0; p < pMax; ++p) {
            prof[p]   = hProf[type][p]->ProfileY(SF("profName%d%d",type, rand()), ii, ii, "i");
            profBB[p] = hProfBB[type][p]->ProfileY(SF("profNameBB%d%d",type, rand()), ii, ii, "i");
            int jj = 7*6 + 1 - i;
            Nprof[p]   = hProf[type][p]->ProfileY(SF("profName%d%d",type, rand()), jj, jj, "i");
            NprofBB[p] = hProfBB[type][p]->ProfileY(SF("profNameBB%d%d",type, rand()), jj, jj, "i");
        }


        prof[0]->SetLineColor(1);
        //prof[0]->Draw("hist e ][ ");
        prof[0]->Draw("axis");


        for(int p = 1; p < pMax; ++p) {
            prof[p]->SetLineColor(p+1);
            prof[p]->Draw("hist e same ][");
        }
        //prof[0]->Draw("hist e same  ][");

        TH1D *myHist;
        if(doSigma) {
            profBB[0]->SetErrorOption("s");
            myHist = new TH1D(SF("myHist%d",rand()), "", profBB[0]->GetNbinsX(), profBB[0]->GetXaxis()->GetXbins()->GetArray());
            for(int k = 1; k <= profBB[0]->GetNbinsX(); ++k) {
                double v = profBB[0]->GetBinContent(k);
                double c = profBB[0]->GetBinCenter(k);
                myHist->SetBinContent(k, (v - c)/c);
                myHist->SetBinError(k, profBB[0]->GetBinError(k)/c);
            }
            
            myHist->SetLineColor(kRed);
            myHist->SetMarkerColor(kRed);
            myHist->Draw("hist e same  ][");
        }
        else {
            profBB[0]->SetLineColor(kRed);
            profBB[0]->SetMarkerColor(kRed);
            profBB[0]->Draw("hist e same  ][");

            NprofBB[0]->SetLineColor(kBlack);
            NprofBB[0]->SetLineWidth(1);
            NprofBB[0]->Draw("hist e same  ][");
        }


        GetXaxis()->SetRangeUser(74, 2000);

        if(doSigma) {
            GetYaxis()->SetRangeUser(-0.15, 0.15);
            if(type == 3 || type == 7)
                GetYaxis()->SetRangeUser(-0.30, 0.30);

            gPad->SetLogy(0);
        }
        else if(type == 3 || type == 7) {
            GetYaxis()->SetRangeUser(74, 2000);
            gPad->SetLogy();
        }
        else {
            if(type == 2 || type == 10)
                GetYaxis()->SetRangeUser(0.91, 1.09);
            else if(type == 4 || type == 5 || type == 6)
                GetYaxis()->SetRangeUser(0.951, 1.09);
            else
                GetYaxis()->SetRangeUser(0.21, 1.25);
            gPad->SetLogy(0);

        }

        GetYaxis()->SetNdivisions(505);
        //GetXaxis()->SetMoreLogLabels();
        //GetYaxis()->SetMoreLogLabels();
        gPad->SetLogx();


        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});
        GetFrame()->SetTitle("");

        double l = hProf[type][0]->GetXaxis()->GetBinLowEdge(ii);
        double u = hProf[type][0]->GetXaxis()->GetBinUpEdge(ii);
        DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u), (type == 3||type ==7) ? 8 : -1);

        if(i == 1) {
            //GetYaxis()->SetTitle("#LT p_{T}^{PUPPI}#GT/#LT p_{T}^{CHS}#GT");

            if(doSigma) {
                if(type == 3) DrawLatexLeft(4.4, "(p_{T}^{PUPPI}(probe) - p_{T}^{CHS})/ p_{T}^{CHS}", -1, "t>");
                else if(type == 7) DrawLatexLeft(4.4, "(p_{T}^{CHS}(probe) - p_{T}^{CHS})/ p_{T}^{CHS}", -1, "t>");

                else if(type == 8) DrawLatexLeft(4.4, "(p_{T}^{CHS} - #bar{p}_{T}^{PUPPI})/ #bar{p}_{T}^{PUPPI}", -1, "t>");
                else if(type == 9) DrawLatexLeft(4.4, "(p_{T}^{PUPPI} - #bar{p}_{T}^{CHS})/ #bar{p}_{T}^{CHS}", -1, "t>");

            }
            else {
                if(type == 3)      DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI} (probe) #GT [GeV]", -1, "t>");
                else if(type == 7) DrawLatexLeft(4.4, "#LT p_{T}^{CHS} (probe) #GT [GeV]", -1, "t>");
                else if(type == 0) DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI} (probe) / p_{T}^{CHS} (tag) #GT", -1, "t>");
                else if(type == 1) DrawLatexLeft(4.4, "#LT p_{T}^{CHS} (probe) / p_{T}^{CHS} (tag) #GT", -1, "t>");
                else if(type == 2 || type == 4 || type == 5 || type == 6 || type == 10)
                                   DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI} (probe) / p_{T}^{CHS} (probe) #GT", -1, "t>");
            }
        }
        if(i == 42) {
            //GetXaxis()->SetTitle("p_{Ttag}^{CHS} [GeV]");
            if(type == 4)
                DrawLatexDown(3.0, "p_{T}^{PUPPI} (probe) [GeV]", -1, "r");
            else if(type == 5)
                DrawLatexDown(3.0, "p_{T}^{CHS} (probe) [GeV]", -1, "r");
            else if(type == 6)
                DrawLatexDown(3.0, "#sqrt{p_{T}^{CHS} p_{T}^{PUPPI}} (probe) [GeV]", -1, "r");

            else if(type == 8) DrawLatexDown(3.0, "p_{T}^{PUPPI} [GeV]", -1, "r");
            else if(type == 9) DrawLatexDown(3.0, "p_{T}^{CHS} [GeV]", -1, "r");
            else if(type == 10) DrawLatexDown(3.0, "p_{T}^{PUPPI} (tag) [GeV]", -1, "r");

            else
                DrawLatexDown(3.0, "p_{T}^{CHS} (tag) [GeV]", -1, "r");

        }

        if( (i == 41 && (type == 3 || type == 7)) || (i == 42 && (type != 3  && type != 7))) {
            TLegend *leg = new TLegend(0.1, 0.55, 0.55, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.8*GetXaxis()->GetTitleSize());
            //leg->AddEntry(prof[0], "any #alpha", "ple");
            if(doSigma) {
                leg->AddEntry(myHist, "#alpha < 0.3 (#eta > 0)", "ple");
            }
            else {
                leg->AddEntry(profBB[0], "#alpha < 0.3 (#eta > 0)", "ple");
                leg->AddEntry(NprofBB[0], "#alpha < 0.3 (#eta < 0)", "ple");
            }
            leg->Draw();
        }


    }

    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}

void PLOTTER::CompareSigmas()
{
    bool doSigma = true;

    TCanvas *can = new TCanvas("can", "can", 800, 500);

    SetLeftRight(0.08, 0.03);

    DividePad(vector<double>(7,1.), vector<double>(6,1.));

    gStyle->SetOptStat(0);

    for(int i = 1; i <= 7*6; ++i) {
        int ii = 7*6 + i;
        can->cd(i);
        TProfile *profBB3;
        TProfile  *NprofBB3;
        TProfile *profBB7;
        TProfile  *NprofBB7;

        int jj = 7*6 + 1 - i;
        int type = 3;
        profBB3 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), ii, ii, "i");
        NprofBB3 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), jj, jj, "i");
        type = 7;
        profBB7 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), ii, ii, "i");
        NprofBB7 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), jj, jj, "i");


        profBB3->SetLineColor(1);
        //prof[0]->Draw("hist e ][ ");
        profBB3->Draw("axis");


        TH1D *myHist;
        if(doSigma) {
            profBB3->SetErrorOption("s");
            profBB7->SetErrorOption("s");
            myHist = new TH1D(SF("myHist%d",rand()), "", profBB3->GetNbinsX(), profBB3->GetXaxis()->GetXbins()->GetArray());
            for(int k = 1; k <= profBB3->GetNbinsX(); ++k) {
                double v3 = profBB3->GetBinContent(k);
                double c3 = profBB3->GetBinCenter(k);
                double r =  profBB7->GetBinError(k) !=  0  ?  profBB3->GetBinError(k) /  profBB7->GetBinError(k)  : 0;
                myHist->SetBinContent(k, r );
                //cout << "HELENKA " << k <<" "<< r << endl;
            }
            
            myHist->SetLineColor(kRed);
            myHist->SetMarkerColor(kRed);
            myHist->Draw("hist e same  ][");
        }


        GetXaxis()->SetRangeUser(74, 2000);

        GetYaxis()->SetRangeUser(0.94, 1.16);
        gPad->SetLogy(0);

        GetYaxis()->SetNdivisions(505);
        //GetXaxis()->SetMoreLogLabels();
        //GetYaxis()->SetMoreLogLabels();
        gPad->SetLogx();


        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});
        GetFrame()->SetTitle("");

        double l = hProf[type][0]->GetXaxis()->GetBinLowEdge(ii);
        double u = hProf[type][0]->GetXaxis()->GetBinUpEdge(ii);
        DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u),  -1);

        if(i == 1) {
            //GetYaxis()->SetTitle("#LT p_{T}^{PUPPI}#GT/#LT p_{T}^{CHS}#GT");
            DrawLatexLeft(4.4, "#sigma^{PUPPI+x}/ #sigma^{CHS+x}", -1, "t>");
        }

        if(i == 42) {
            DrawLatexDown(3.0, "p_{T}^{CHS} (tag) [GeV]", -1, "r");

        }


        if(i == 42) {
            TLegend *leg = new TLegend(0.1, 0.55, 0.55, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.8*GetXaxis()->GetTitleSize());
            //leg->AddEntry(prof[0], "any #alpha", "ple");
            leg->AddEntry(myHist, "#alpha < 0.3 (pos #eta)", "ple");
            leg->Draw();
        }



    }

    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}

void PLOTTER::PlotVar1overVar2()
{
    TCanvas *can = new TCanvas("can", "can", 800, 500);

    SetLeftRight(0.08, 0.03);

    DividePad(vector<double>(7,1.), vector<double>(6,1.));

    gStyle->SetOptStat(0);

    for(int i = 1; i <= 7*6; ++i) {
        int ii = 7*6 + i;
        can->cd(i);
        TProfile *profBB1;
        TProfile  *NprofBB1;
        TProfile *profBB2;
        TProfile  *NprofBB2;

        int jj = 7*6 + 1 - i;
        int type = 0;
        profBB1 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), ii, ii, "i");
        NprofBB1 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), jj, jj, "i");
        type = 1;
        profBB2 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), ii, ii, "i");
        NprofBB2 = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), jj, jj, "i");


        profBB1->SetLineColor(1);
        //prof[0]->Draw("hist e ][ ");
        profBB1->Draw("axis");


        //TH1D *profP = dynamic_cast<TH1D*>(profBB1->Clone(SF("%d", rand())));
        //TH1D *profN = dynamic_cast<TH1D*>(NprofBB1->Clone(SF("%d",rand())));

        TH1D *profP = new TH1D(SF("myHistP%d",rand()), "", profBB1->GetNbinsX(), profBB1->GetXaxis()->GetXbins()->GetArray());
        TH1D *profN = new TH1D(SF("myHistN%d",rand()), "", profBB1->GetNbinsX(), profBB1->GetXaxis()->GetXbins()->GetArray());


        for(int k = 1; k < profP->GetNbinsX(); ++k) {
            double rP = profBB2->GetBinContent(k) != 0 ?  profBB1->GetBinContent(k) /  profBB2->GetBinContent(k) : 0;
            double rN = NprofBB2->GetBinContent(k) != 0 ?  NprofBB1->GetBinContent(k) /  NprofBB2->GetBinContent(k) : 0;
            cout << "RADECEK " << rP << endl;
            profP->SetBinContent(k, rP);
            profP->SetBinError(k, 0);
            profN->SetBinContent(k, rN);
            profN->SetBinError(k, 0);
        }

        profP->SetLineColor(kRed);
        profP->SetMarkerColor(kRed);
        profP->Draw("hist e same  ][");

        profN->SetLineColor(kBlack);
        profN->SetMarkerColor(kBlack);
        profN->Draw("hist e same  ][");

        GetXaxis()->SetRangeUser(74, 2000);

        //GetYaxis()->SetRangeUser(0.94, 1.16);
        GetYaxis()->SetRangeUser(0.91, 1.09);
        gPad->SetLogy(0);

        GetYaxis()->SetNdivisions(505);
        //GetXaxis()->SetMoreLogLabels();
        //GetYaxis()->SetMoreLogLabels();
        gPad->SetLogx();


        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});
        GetFrame()->SetTitle("");

        double l = hProf[type][0]->GetXaxis()->GetBinLowEdge(ii);
        double u = hProf[type][0]->GetXaxis()->GetBinUpEdge(ii);
        DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u),  -1);


        if(i == 1) {
            //GetYaxis()->SetTitle("#LT p_{T}^{PUPPI}#GT/#LT p_{T}^{CHS}#GT");
            DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI} (probe) / p_{T}^{CHS} (tag) #GT  /  #LT p_{T}^{CHS} (probe) / p_{T}^{CHS} (tag) #GT", -1, "t>");

                //else if(type == 0) DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI} (probe) / p_{T}^{CHS} (tag) #GT", -1, "t>");
                //else if(type == 1) DrawLatexLeft(4.4, "#LT p_{T}^{CHS} (probe) / p_{T}^{CHS} (tag) #GT", -1, "t>");

        }

        if(i == 42) {
            DrawLatexDown(3.0, "p_{T}^{CHS} (tag) [GeV]", -1, "r");
        }


        if(i == 42) {
            TLegend *leg = new TLegend(0.1, 0.55, 0.55, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.8*GetXaxis()->GetTitleSize());
            //leg->AddEntry(prof[0], "any #alpha", "ple");
            leg->AddEntry(profP, "#alpha < 0.3 (#eta > 0)", "ple");
            leg->AddEntry(profN, "#alpha < 0.3 (#eta < 0)", "ple");
            leg->Draw();
        }
        //*/



    }

    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}

















void PLOTTER::PlotProfilesOverlap()
{
    TCanvas *can = new TCanvas("can", "can", 800, 500);
    SetLeftRight(0.08, 0.03);
    DividePad(vector<double>(7,1.), vector<double>(6,1.));
    gStyle->SetOptStat(0);

    for(int i = 1; i <= 7*6; ++i) {
        can->cd(i);
        TProfile *prof[8], *profBB[8];
        for(int type = 0; type < 7; ++type) {
            prof[type]   = hProf[type][0]->ProfileY(SF("profName%d%d",type, rand()), i, i, "i");
            profBB[type] = hProfBB[type][0]->ProfileY(SF("profNameBB%d%d",type, rand()), i, i, "i");
        }

        auto plot = [&](int i, int c, TString tag) {
            prof[i]->SetLineColor(c);
            prof[i]->SetMarkerColor(c);
            prof[i]->Draw(tag);
        };
        plot(2,1, "e");
        plot(4,2, "e same");
        plot(5,3, "e same");

        GetXaxis()->SetRangeUser(74, 2000);
        GetYaxis()->SetRangeUser(0.951, 1.09);

        GetYaxis()->SetNdivisions(505);
        //GetXaxis()->SetMoreLogLabels();
        gPad->SetLogx();


        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});
        GetFrame()->SetTitle("");

        double l = hProf[2][0]->GetXaxis()->GetBinLowEdge(i);
        double u = hProf[2][0]->GetXaxis()->GetBinUpEdge(i);
        DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u));

        if(i == 1) 
           DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI} (probe) / p_{T}^{CHS} (probe) #GT", -1, "t>");
        if(i == 42)
            DrawLatexDown(3.0, "p_{T} [GeV]", -1, "r");


        if(i == 42) {
            TLegend *leg = new TLegend(0.17, 0.45, 0.65, 0.85);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.8*GetXaxis()->GetTitleSize());
            leg->AddEntry(prof[2], "p_{T}^{CHS} (tag)", "ple");
            leg->AddEntry(prof[4], "p_{T}^{PUPPI} (probe)", "ple");
            leg->AddEntry(prof[5], "p_{T}^{CHS} (probe)", "ple");
            leg->Draw();
        }


    }

    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}






pair<TGraphErrors*, TGraphErrors*> GetMikkoCorr(vector<TProfile2D*> *hProf, int i)
{
    TProfile  *prof[4];
    for(int type = 0; type < 4; ++type)
        prof[type] = hProf[type][0]->ProfileY(SF("profName%d%d",type, rand()), i, i, "i");

    TGraphErrors *gr1 = new TGraphErrors();
    TGraphErrors *gr2 = new TGraphErrors();

    for(int k = 1; k <= prof[0]->GetNbinsX(); ++k) {
        double v1 = prof[0]->GetBinContent(k);
        double v2 = prof[1]->GetBinContent(k);
        double v3 = prof[2]->GetBinContent(k);
        double e3 = prof[2]->GetBinError(k);
        double v4 = prof[3]->GetBinContent(k);
        double e4 = prof[3]->GetBinError(k);
        double r = v2 != 0 ? v1 / v2 : 0;
        gr1->SetPoint(k-1, v4, r); 
        gr1->SetPointError(k-1, e4, e3); 

        gr2->SetPoint(k-1, v4, v3); 
        gr2->SetPointError(k-1, e4, e3); 
        //cout << k << " " << v4 << " "<< r << endl;
    }
    return make_pair(gr1, gr2);

}

void CleanGraph(TGraphErrors *gr)
{
    double xOrg = -30000;
    int iLast = 1000;
    for(int i = 0; i < gr->GetN(); ++i) {
        double x, y;
        gr->GetPoint(i, x, y);
        if(x < xOrg) {
            iLast = i;
            break;
        }
        xOrg = x;
    }
    iLast = max(0, iLast);
    //cout << "Old Npoints " << gr->GetN() << endl;

    for(int i = gr->GetN()-1; i >= iLast; --i) {
        //cout << "Cleaning " << i << endl;
        gr->RemovePoint(i);
    }
    //cout << "New Npoints " << gr->GetN() << endl;
}

pair<TH1D*,TH1D*> PLOTTER::PlottMatchingCorr(int direction, bool doFit = true)
{
    TCanvas *can = new TCanvas("can", "can", 800, 500);
    SetLeftRight(0.08, 0.03);

    DividePad(vector<double>(7,1.), vector<double>(6,1.));


    TH1D *hFit84 = nullptr, *hFit1000 = nullptr;
    if(doFit) {
        const double *ar = hBalEtaPtAll[0][0]->GetXaxis()->GetXbins()->GetArray();
        vector<double> arr;
        for(int k = 0; k <= hBalEtaPtAll[0][0]->GetNbinsX(); ++k) {
            cout << "Holka " << ar[k] << endl;
            if(ar[k] >= -0.00001) arr.push_back(ar[k]);
        }

        hFit84   =  new TH1D("hFit84", "fit84",  arr.size()-1, arr.data());
        hFit1000 =  new TH1D("hFit1000", "fit1000",  arr.size()-1, arr.data());
    }



    gStyle->SetOptStat(0);

    for(int i = 1; i <= 7*6; ++i) {
        int ii = 7*6+i;
        if(direction == -1)
            ii = 7*6+1-i;

        can->cd(i);

        cout << "Helenka " << __LINE__ << endl;
        TGraphErrors *gr1, *gr2;
        tie(gr1, gr2) = GetMikkoCorr(hProf, ii);
        TGraphErrors *gr1BB, *gr2BB;
        tie(gr1BB, gr2BB) = GetMikkoCorr(hProfBB, ii);

        hBalEtaPtAll[0][0]->GetXaxis()->SetRange(ii,ii);
        TH2D *hTemp = dynamic_cast<TH2D*>(hBalEtaPtAll[0][0]->Project3D(SF("%d_yz",rand()))); 
        assert(hTemp);
        TProfile *profOrg = hTemp->ProfileY();

        profOrg->Draw("hist");

        cout << "Helenka " << __LINE__ << endl;
        CleanGraph(gr1BB);
        //gr1->Draw("pe");
        gr1BB->SetLineColor(kBlack);
        gr1BB->SetMarkerColor(kBlack);
        gr1BB->Draw("pe same");

        GetXaxis()->SetRangeUser(74, 2000);
        //GetXaxis()->SetMoreLogLabels();
        GetYaxis()->SetMoreLogLabels();
        GetYaxis()->SetRangeUser(0.951, 1.19);
        gPad->SetLogx();

        cout << "Helenka " << __LINE__ << endl;
        GetYaxis()->SetNdivisions(505);

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});

        GetFrame()->SetTitle("");
        if(i == 1) {
            //GetYaxis()->SetTitle("#LT p_{T}^{PUPPI}#GT/#LT p_{T}^{CHS}#GT");
            DrawLatexLeft(4.4, "#LT p_{T}^{PUPPI}#GT/#LT p_{T}^{CHS}#GT", -1, "t>");
        }
        if(i == 42) {
            GetXaxis()->SetTitle("p_{T}^{PUPPI} [GeV]");

            TLegend *leg = new TLegend(0.1, 0.5, 0.65, 0.85);
            leg->SetBorderSize(0);
            //leg->AddEntry(gr1, "any #alpha", "lpe");
            leg->AddEntry(gr1BB, "#alpha < 0.3", "lpe");
            leg->AddEntry(profOrg, "Old Method", "l");
            leg->Draw();
        }

        cout << "Helenka " << __LINE__ << endl;


        double l = hBalEtaPtAll[0][0]->GetXaxis()->GetBinLowEdge(ii);
        double u = hBalEtaPtAll[0][0]->GetXaxis()->GetBinUpEdge(ii);
        double c = hBalEtaPtAll[0][0]->GetXaxis()->GetBinCenter(ii);

        cout << "Helenka " << __LINE__ << endl;

        /////////////////////////////////////////////////////////////
        double y = hBalEtaPtAll[0][0]->GetXaxis()->GetBinCenter(ii);
        cout << "my y " << y << endl;
        TF1 *fun = new TF1(SF("fun%d",rand()), "[0] + [1]*log(x/84)", 84, 2000);
        fun->SetParameter(0, 1);
        fun->SetParameter(1, 0);
        if(doFit && 1) {
            double maxPt = 0;
            cout << "Helenka " << __LINE__ << endl;
            for(int i = 0; i <= gr1BB->GetN(); ++i) {
                double x, y;
                gr1BB->GetPoint(i, x, y);
                maxPt = max(maxPt, x);
            }
            double maxPtNow = std::min(500.0, maxPt);
            if(maxPtNow > 84) {
                cout << "Helenka " << __LINE__ << endl;

                gr1BB->Fit(fun, "", "", 84, maxPtNow);
                cout << "HuHu " << 84 << " "<< maxPtNow << endl;
                fun->DrawF1(84, maxPtNow, "l same");
                TF1 *fun2 = (TF1*) fun->Clone(SF("fun2_%d",i));
                fun2->SetLineColor(kGreen);
                if(maxPt > maxPtNow)
                    fun2->DrawF1(maxPtNow, maxPt,  "l same");

                cout << "Helenka " << __LINE__ << endl;


                double yMin = hBalEtaPtAll[0][0]->GetXaxis()->GetBinLowEdge(ii);
                double yMax = hBalEtaPtAll[0][0]->GetXaxis()->GetBinUpEdge(ii);
                cout << "Radek "<<  yMin << " " << yMax <<
                    " 84 "<< maxPt<<" 84 "<< maxPtNow<<" "<< /*round(6500/cosh(yMin))<<" "<<*/ round(100000*fun->GetParameter(0))/100000. <<" "<< round(100000*fun->GetParameter(1))/100000. << endl;


                double p0 = fun->GetParameter(0);
                double p1 = fun->GetParameter(1);
                double p0err = fun->GetParError(0);
                double p1err = fun->GetParError(1);

                double s84=0, ss84=0;
                double s1000=0, ss1000=0;
                int N = 10000;
                for(int k = 0; k < N; ++k) {
                    double p0n = p0 + gRandom->Gaus(p0, p0err);
                    double p1n = p1 + gRandom->Gaus(p1, p1err);
                    fun->SetParameters(p0n, p1n);
                    s84  += fun->Eval(84);
                    ss84 += pow(fun->Eval(84),2);
                    s1000  += fun->Eval(1000);
                    ss1000 += pow(fun->Eval(1000),2);
                }
                fun->SetParameters(p0, p1);
                cout << "Helenka " << __LINE__ << endl;

                //cout << "HEL84 " << ss84/ << 
                double e84   = sqrt(ss84/N - pow(s84/N,2));
                double e1000 = sqrt(ss1000/N - pow(s1000/N,2));
                cout <<"HELAX "<< e84 <<" "<< e1000 << endl;

                cout << "Helenka " << __LINE__ << endl;

                int iBin = hFit84->FindBin(abs(c));
                if(maxPtNow > 84 && abs(c) < 4.4) {
                    hFit84->SetBinContent(iBin, fun->Eval(84));
                    hFit84->SetBinError(iBin, e84);
                    hFit1000->SetBinContent(iBin, fun->Eval(1000));
                    hFit1000->SetBinError(iBin, e1000);
                }

            }

            cout << "Helenka " << __LINE__ << endl;

        }


        cout << "Helenka " << __LINE__ <<" "<< l <<" "<< u<< endl;
        DrawLatexUp(-1, SF("%1.2f < #eta < %1.2f", l, u));

        cout << "Helenka " << __LINE__ << endl;


    }

        cout << "Helenka " << __LINE__ << endl;
    can->Print(outName);
        cout << "Helenka " << __LINE__ << endl;
    can->Clear();
    delete can;
    return make_pair(hFit84, hFit1000);
}

void  PLOTTER::PlotFit(TH1D *h84, TH1D *h1000)
{
    TCanvas *can = new TCanvas("can", "can", 500, 500);
    SetLeftRight(0.15, 0.15);

    h84->SetLineColor(kBlack);
    h84->Draw("][");
    h1000->SetLineColor(kBlue);
    h1000->Draw("][ same");

    GetYaxis()->SetRangeUser(0.92, 1.08);
    GetXaxis()->SetRangeUser(0.0, 4.35);
    GetYaxis()->SetTitle("p_{T}^{PUPPI}/p_{T}^{CHS} Correction");
    GetXaxis()->SetTitle("|#eta|");
    GetFrame()->SetTitle("");

    TLegend *leg = new TLegend(0.2, 0.15, 0.5, 0.4);
    leg->SetBorderSize(0);
    leg->SetTextSize(GetXaxis()->GetTitleSize());
    leg->AddEntry((TObject*)nullptr, "The PUPPI/CHS", "h");
    leg->AddEntry((TObject*)nullptr, "matching corrections", "h");
    leg->AddEntry((TObject*)nullptr, "(logLin fit)", "h");
    leg->AddEntry(h84, "p_{T}^{PUPPI} = 84 GeV", "l");
    leg->AddEntry(h1000, "p_{T}^{PUPPI} = 1000 GeV", "l");
    leg->Draw();




    DrawLatexUp(1, "2016 RunB-H (13 TeV)", -1, "vr" );

    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}

void  PLOTTER::PlotFitFour(TH1D *h84p, TH1D *h1000p, TH1D *h84n, TH1D *h1000n)
{
    TCanvas *can = new TCanvas("can", "can", 500, 500);
    SetLeftRight(0.15, 0.15);

    h84p->SetLineColor(kRed);
    h84p->Draw("hist  e ][");
    h84p->Draw("hist   ][ same");
    h1000p->SetLineColor(kRed);
    h1000p->SetLineStyle(2);
    h1000p->Draw("hist e ][ same");
    h1000p->Draw("hist   ][ same");

    h84n->SetLineColor(kBlack);
    h84n->Draw("hist e ][ same");
    h84n->Draw("hist   ][ same");
    h1000n->SetLineColor(kBlack);
    h1000n->SetLineStyle(2);
    h1000n->Draw("hist e ][ same");
    h1000n->Draw("hist   ][ same");


    GetYaxis()->SetRangeUser(0.92, 1.08);
    GetXaxis()->SetRangeUser(0.0, 4.35);
    GetYaxis()->SetTitle("p_{T}^{PUPPI}/p_{T}^{CHS} Correction");
    GetXaxis()->SetTitle("|#eta|");
    GetFrame()->SetTitle("");

    TLegend *leg = new TLegend(0.2, 0.15, 0.5, 0.4);
    leg->SetBorderSize(0);
    leg->SetTextSize(GetXaxis()->GetTitleSize());
    leg->AddEntry((TObject*)nullptr, "The PUPPI/CHS", "h");
    leg->AddEntry((TObject*)nullptr, "matching corrections", "h");
    leg->AddEntry((TObject*)nullptr, "(logLin fit)", "h");
    //leg->AddEntry(h84p, "p_{T}^{PUPPI} = 84 GeV", "l");
    //leg->AddEntry(h1000p, "p_{T}^{PUPPI} = 1000 GeV", "l");
    leg->Draw();


    TLegend *leg2 = new TLegend(0.2, 0.65, 0.7, 0.8);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.7*GetXaxis()->GetTitleSize());
    leg2->SetNColumns(2);
    leg2->AddEntry((TObject*)nullptr, "pos", "");
    leg2->AddEntry((TObject*)nullptr, "neg", "");

    leg2->AddEntry(h84p, "p_{T}^{PUPPI} = 84 GeV", "l");
    leg2->AddEntry(h84n, "p_{T}^{PUPPI} = 84 GeV", "l");
    leg2->AddEntry(h1000p, "p_{T}^{PUPPI} = 1000 GeV", "l");
    leg2->AddEntry(h1000n, "p_{T}^{PUPPI} = 1000 GeV", "l");
    leg2->Draw();



    DrawLatexUp(1, "2016 RunB-H (13 TeV)", -1, "vr" );

    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}







pair<TPad*, TPad*> TitleSpace(int perID)
{
    TVirtualPad *can = gPad;
    const double sp = 0.98;
    TPad *up = new TPad(SF("%d",rand()),"",   0, sp, 1, 1.0);

    up->Draw();

    up->cd();
    TLatex *lat = new TLatex;
    TString s;
    if(perID != 0) s =  SF("Run 2016%c", perID + 'A');
    else s =  SF("Run 2016 all");
    lat->SetTextSize(1.0);
    lat->DrawLatexNDC(0.8, 0.0, s);


    can->cd();


    TPad *dn = new TPad(SF("%d",rand()),"", 0, 0.0, 1, sp);
    dn->Draw();
    dn->cd();
    //return make_pair(up, dn);
    return make_pair(up, dn);
}

void SetFontSizes(TH1 *h, double val)
{
    h->GetXaxis()->SetLabelSize(val);
    h->GetYaxis()->SetLabelSize(val);
    h->GetXaxis()->SetTitleSize(val);
    h->GetYaxis()->SetTitleSize(val);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
}




void PLOTTER::JEC()
{
    TCanvas *can = new TCanvas("can", "can");

    //hJECpuppi[perID]->
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    auto hJECp = hJECpuppi[perID];
    auto hJECc = hJECchs[perID];

    const vector<int> etaBins = {1, 4, 7, 10, 13, 16};
    const vector<int> ptBins = {6, 12, 18, 24, 30, 36};

    dnPad->Divide(5, 5, 0, 0);
    dnPad->SetBottomMargin(0.1);
    for(int i = 1; i <= 5; ++i) 
    for(int j = 1; j <= 5; ++j) {
        dnPad->cd((i-1)*5 + j);
        int iS = etaBins[i-1];
        int iE = etaBins[i+1-1]-1;
        int jS = ptBins[j];
        int jE = ptBins[j+1]-1;

        TH1D *hPuppi  = hJECp->ProjectionZ(SF("ProjBoth%d%d%d",i,j,rand()),  iS, iE, jS, jE);
        TH1D *hChs    = hJECc->ProjectionZ(SF("ProjBothUp%d%d%d",i,j,rand()),  iS, iE, jS, jE);

        //cout <<"Radek "<< i << " "<< j<<" " << hBoth->GetStdDev() << endl;

        double lPt = hJECp->GetYaxis()->GetBinLowEdge(jS);
        double uPt = hJECp->GetYaxis()->GetBinUpEdge(jE);
        double lEt = hJECp->GetXaxis()->GetBinLowEdge(iS);
        double uEt = hJECp->GetXaxis()->GetBinUpEdge(iE);

        //cout << "RADEK " << j <<" "<< lPt << " "<< endl;

        //hBoth->SetTitle(SF("#splitline{%3.0f < p_{T} < %3.0f}{%1.1f < #eta < %1.1f}", lPt, uPt, lEt, uEt));
        //cout << "RADEKT " << SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt) << endl;
        hChs->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));

        hChs->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hChs->SetMaximum(1.2*hChs->GetMaximum());
        hChs->SetLineColor(kBlack);
        SetFontSizes(hChs, 0.07);
        hChs->SetTitleSize(0.13);
        hChs->Draw("hist e ");

        hPuppi->SetLineColor(kBlue);
        hPuppi->SetLineStyle(1);


        hPuppi->Draw("hist e same");

        TLegend *leg = new TLegend(0.8, 0.5, 0.95, 0.65);
        leg->SetBorderSize(0);
        leg->AddEntry(hChs, "CHS");
        leg->AddEntry(hPuppi, "PUPPI");
        leg->Draw();
        //hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));
    }



    can->Print(outName);
    can->Clear();
    delete can;
}


void PLOTTER::PtEtaDep()
{
    auto hPUPPI = hEtaPtPUPPI[perID];
    auto hPUPPIa = hEtaPtPUPPIalone[perID];
    auto hCHS = hEtaPtCHS[perID];

    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //can->
    //dnPad->Divide(5, 4, 0.0001, 0.002);
    dnPad->cd();
    DividePad(vector<double>(5,1.), vector<double>(4,1.));
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 20; ++i) {
        dnPad->cd(i);
        //gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
        gPad->SetLogx();
        //gPad->SetLogy();
        //if(i >= 19) continue;
        int iS = 1*i ;
        int iE = 1*(i+1) -1 ;
        TH1D *hRef   = hPUPPI->ProjectionY(SF("ProjEta%d%d",i,rand()),  iS, iE);
        TH1D *hPuppi = hPUPPI->ProjectionY(SF("ProjEta%d%d",i,rand()),  iS, iE);
        TH1D *hPuppiA = hPUPPIa->ProjectionY(SF("ProjEtaAll%d%d",i,rand()),  iS, iE);
        TH1D *hChs    = hCHS->ProjectionY(SF("ProjEtaCHS%d%d",i,rand()),  iS, iE);
        for(int i = 0; i <= hRef->GetNbinsX(); ++i)
            hRef->SetBinError(i,0);

        double l = hPUPPI->GetXaxis()->GetBinLowEdge(iS);
        double u = hPUPPI->GetXaxis()->GetBinUpEdge(iE);
        //cout <<"Helenka " << l <<" "<< u << endl;
        hPuppiA->SetLineColor(kBlack);
        hChs->SetLineColor(kRed);


        //hPuppi->SetMaximum(1.2*hPUPPI->GetMaximum());
        //hPuppi->SetMaximum(10);

        //hPuppi->Scale(1., "width");
        hPuppiA->Divide(hRef);
        hPuppi->Divide(hRef);
        hChs->Divide(hRef);

        //hChs->Scale(1., "width");

        //hPuppiA->SetMinimum(0.00);
        //hPuppiA->SetMaximum(0.20);

        hPuppiA->Draw();
        hPuppi->Draw("same");
        hChs->Draw("same");
        //SetFontSizes(hPuppiA, 0.07);

        GetXaxis()->SetRangeUser(34, 1780);
        GetYaxis()->SetRangeUser(0.41, 1.59);
        GetXaxis()->SetMoreLogLabels();
        GetXaxis()->SetNoExponent();

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 3.5});

        GetFrame()->SetTitle("");
        if(i == 1)
            GetYaxis()->SetTitle("#sigma/#sigma_{PUPPI}");
        if(i == 20)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} [GeV]");

        DrawLatexUp(-1, SF("%1.1f < |#eta| < %1.1f", l, u));

        if(i == 19) {
            TLegend *leg = new TLegend(0.2, 0.5, 0.7, 0.7);
            leg->SetBorderSize(0);
            leg->SetTextSize(GetXaxis()->GetTitleSize());
            leg->AddEntry(hPuppi, "PUPPI");
            leg->AddEntry(hChs, "CHS");
            leg->Draw();
        }

    }
    
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}

void PLOTTER::Unmatched(TString type = "tested")
{
    TH2D* hPUPPI = nullptr;
    TH2D* hPUPPIa = nullptr;

    if(type == "tested") {
        hPUPPI = hEtaPtPUPPI[perID];
        hPUPPIa = hEtaPtPUPPIalone[perID];
    }
    else {
        hPUPPI = hEtaPtCHS[perID];
        hPUPPIa = hEtaPtCHSalone[perID];
    }

    TCanvas *can = new TCanvas("can", "can", 800, 500);

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //can->
    //dnPad->Divide(5, 4, 0.0001, 0.002);
    dnPad->cd();
    SetLeftRight(0.10, 0.03);
    DividePad(vector<double>(7,1.), vector<double>(6,1.));
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 42; ++i) {
        dnPad->cd(i);
        //gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
        gPad->SetLogx();
        //gPad->SetLogy();
        //if(i >= 19) continue;
        int iS = 1*i ;
        int iE = 1*(i+1) -1 ;
        TH1D *hRef   = hPUPPI->ProjectionY(SF("ProjEta%d%d",i,rand()),  iS, iE);
        TH1D *hPuppi = hPUPPI->ProjectionY(SF("ProjEta%d%d",i,rand()),  iS, iE);
        TH1D *hPuppiA = hPUPPIa->ProjectionY(SF("ProjEtaAll%d%d",i,rand()),  iS, iE);
        for(int i = 0; i <= hRef->GetNbinsX(); ++i)
            hRef->SetBinError(i,0);

        double l = hPUPPI->GetXaxis()->GetBinLowEdge(iS);
        double u = hPUPPI->GetXaxis()->GetBinUpEdge(iE);
        //cout <<"Helenka " << l <<" "<< u << endl;
        hPuppiA->SetLineColor(kBlack);

        //hPuppi->SetMaximum(1.2*hPUPPI->GetMaximum());
        //hPuppi->SetMaximum(10);

        //hPuppi->Scale(1., "width");
        hPuppiA->Divide(hPuppi);

        hPuppiA->Print("all");

        //hChs->Scale(1., "width");

        //hPuppiA->SetMinimum(0.00);
        //hPuppiA->SetMaximum(0.20);

        hPuppiA->Draw();
        //hPuppi->Draw("same");
        //hChs->Draw("same");

        //SetFontSizes(hPuppiA, 0.07);

        GetXaxis()->SetRangeUser(34, 1780);
        GetYaxis()->SetRangeUser(0.0, 0.129);
        GetYaxis()->SetNdivisions(505);

        GetXaxis()->SetMoreLogLabels();
        GetXaxis()->SetNoExponent();

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 5.4});

        GetFrame()->SetTitle("");
        if(i == 1)
            GetYaxis()->SetTitle("#sigma^{unmatched}_{PUPI}/#sigma_{PUPPI}");
        if(i == 42)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} [GeV]");

        DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u));
    }
    
    //exit(0);
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}


void PLOTTER::MatchingFactorsTimeDep(int rho)
{

    TH3D *hBalEtaPt;
    hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAll[rho][0]->Clone(SF("%d",rand())));
    assert(hBalEtaPt);


    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(0);

    //can->
    //dnPad->Divide(5, 4, 0.0001, 0.002);
    dnPad->cd();
    DividePad(vector<double>(5,1.), vector<double>(4,1.));
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 20; ++i) {
        dnPad->cd(i);
        //gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
        gPad->SetLogx();
        //gPad->SetLogy();
        //if(i >= 19) continue;
        int iS = 1*i ;
        int iE = 1*(i+1) -1 ;


        double l = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);

        //cout << "xOrg " << hBalEtaPt->GetXaxis()->GetXmin() <<" "<< hBalEtaPt->GetXaxis()->GetXmax()<< endl;
        //cout << "yOrg " << hBalEtaPt->GetYaxis()->GetXmin() <<" "<< hBalEtaPt->GetYaxis()->GetXmax()<< endl;
        //cout << "zOrg " << hBalEtaPt->GetZaxis()->GetXmin() <<" "<< hBalEtaPt->GetZaxis()->GetXmax()<< endl;

        vector<TProfile*> prof(8);
        for(int k = 0; k < 8; ++k) {
            iS = min(19,iS);
            iE = min(19,iE);
            hBalEtaPtAll[rho][k]->GetXaxis()->SetRange(iS,iE);
            TH2D *hTemp = dynamic_cast<TH2D*>(hBalEtaPtAll[rho][k]->Project3D(SF("%d_yz",rand()))); 
            assert(hTemp);
            prof[k] = hTemp->ProfileY();

            prof[k]->SetLineColor(1+k);
            if(k==0)
                prof[k]->Draw("hist");
            else
                prof[k]->Draw("hist same");
        }
        prof[0]->Draw("hist same");

        GetXaxis()->SetRangeUser(34, 1780);
        GetYaxis()->SetRangeUser(0.901, 1.149);
        GetXaxis()->SetMoreLogLabels();
        GetXaxis()->SetNoExponent();
        GetYaxis()->SetNdivisions(505);

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});

        GetFrame()->SetTitle("");
        if(i == 1)
            GetYaxis()->SetTitle("#LT p_{T}^{PUPI}/p_{T}^{CHS}#GT");
        if(i == 20)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} [GeV]");

        DrawLatexUp(-1, SF("%1.1f < |#eta| < %1.1f", l, u));

        if(i == 20) {
            TLegend * leg = new TLegend(0.05, 0.40, 0.32, 0.85);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            for(int k = 0; k < 8; ++k) {
                TString t = k != 0 ? SF("Run %c",'A'+ k) : "All Runs";
                leg->AddEntry(prof[k], t, "l");
            }
            leg->Draw();
        }



    }
    
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}


void CleanProfile(TProfile *prof, double cut)
{
    //Cleaning
    for(int aa = 1; aa <= prof->GetNbinsX(); ++aa) {
        if(prof->GetBinEffectiveEntries(aa) < cut)
            prof->SetBinEntries(aa, 0);
    }
}

void CleanProfile(TProfile2D *prof, double cut)
{
    //Cleaning
    for(int aa = 1; aa <= prof->GetNbinsX() * prof->GetNbinsY(); ++aa) {
        if(prof->GetBinEffectiveEntries(aa) < cut)
            prof->SetBinEntries(aa, 0);
    }
}



pair<TH1D*,TH1D*> PLOTTER::MatchingFactorsPileUpDep(bool plotFit, TString plotType)
{
    int nRuns = 1;
    if(plotType == "time") nRuns = 8;
    else if(plotType == "pileUp") nRuns=5;




    TH3D *hBalEtaPt;
    hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAll[0][0]->Clone(SF("%d",rand())));
    assert(hBalEtaPt);

    TH1D *hFit84 = nullptr, *hFit1000 = nullptr;
    if(plotFit) {
        hFit84   =  new TH1D("hFit84", "fit84",  hBalEtaPt->GetNbinsX(), hBalEtaPt->GetXaxis()->GetXbins()->GetArray());
        hFit1000 =  new TH1D("hFit1000", "fit1000",  hBalEtaPt->GetNbinsX(), hBalEtaPt->GetXaxis()->GetXbins()->GetArray());
    }


    TCanvas *can = new TCanvas("can", "can", 800, 500);

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(0);

    //can->
    //dnPad->Divide(5, 4, 0.0001, 0.002);
    dnPad->cd();
    SetLeftRight(0.08, 0.03);
    DividePad(vector<double>(7,1.), vector<double>(6,1.));
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 7*6; ++i) {
        dnPad->cd(i);
        //gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
        gPad->SetLogx();
        //gPad->SetLogy();
        //if(i >= 19) continue;
        int iS = 1*i ;
        int iE = 1*(i+1) -1 ;


        double l = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);

        //cout << "xOrg " << hBalEtaPt->GetXaxis()->GetXmin() <<" "<< hBalEtaPt->GetXaxis()->GetXmax()<< endl;
        //cout << "yOrg " << hBalEtaPt->GetYaxis()->GetXmin() <<" "<< hBalEtaPt->GetYaxis()->GetXmax()<< endl;
        //cout << "zOrg " << hBalEtaPt->GetZaxis()->GetXmin() <<" "<< hBalEtaPt->GetZaxis()->GetXmax()<< endl;

        vector<TProfile*> prof(8);
        for(int k = 0; k < nRuns; ++k) {
            //iS = min(19,iS);
            //iE = min(19,iE);
            TH3D *ptrNow = (plotType=="time") ? hBalEtaPtAll[0][k] : hBalEtaPtAll[k][0];
            ptrNow->GetXaxis()->SetRange(iS,iE);
            TH2D *hTemp = dynamic_cast<TH2D*>(ptrNow->Project3D(SF("%d_yz",rand()))); 
            assert(hTemp);
            prof[k] = hTemp->ProfileY();

            if(nRuns == 1)
                CleanProfile(prof[k], 100);

            prof[k]->SetLineColor(1+k);
            if(k==0)
                prof[k]->Draw("hist ][");
            else
                prof[k]->Draw("hist same ][");
        }
        prof[0]->Draw("hist e same ][");


        double y = hBalEtaPt->GetXaxis()->GetBinCenter(iS);
        cout << "my y " << y << endl;
        TF1 *fun = new TF1(SF("fun%d",i), "[0] + [1]*log(x/84)", 84, 2000);
        if(plotFit) {

            double maxPt = 0;
            for(int i = prof[0]->FindBin(86); i <= prof[0]->GetNbinsX(); ++i)
                if(prof[0]->GetBinContent(i) < 0.3) {
                    maxPt = prof[0]->GetXaxis()->GetBinLowEdge(i);
                    break;
                }
            double maxPtNow = std::min(500.0, maxPt);

            prof[0]->Fit(fun, "", "", 84, maxPtNow);
            fun->DrawF1(84, maxPtNow, "l same");
            TF1 *fun2 = (TF1*) fun->Clone(SF("fun2_%d",i));
            fun2->SetLineColor(kBlue);
            if(maxPt > maxPtNow)
                fun2->DrawF1(maxPtNow, maxPt,  "l same");


            double yMin = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
            double yMax = hBalEtaPt->GetXaxis()->GetBinUpEdge(iS);
            cout << "Radek "<<  yMin << " " << yMax <<
                    " 84 "<< maxPtNow <<" "<< maxPt <<" "<< round(6500/cosh(yMin))<<" "<< round(1000*fun->GetParameter(0))/1000. <<" "<< round(1000*fun->GetParameter(1))/1000. << endl;

            int iBin = hFit84->FindBin(y);
            if(maxPtNow > 84 && y < 4.4) {
                hFit84->SetBinContent(iBin, fun->Eval(84));
                hFit1000->SetBinContent(iBin, fun->Eval(1000));
            }

        }



        /*
        //Plot result of the fit
        TH1D *hfit1 = dynamic_cast<TH1D*>(hFit->ProjectionX(SF("%d",rand()), iS, iE));
        cout << "I am here " << __LINE__ << prof[0] << endl;
        hfit1->SetLineColor(kRed);
        cout << "I am here " << __LINE__ << prof[0] << endl;
        hfit1->Draw("hist same");
        */



        GetXaxis()->SetRangeUser(84, 2000);
        GetYaxis()->SetRangeUser(0.961, 1.149);
        GetXaxis()->SetMoreLogLabels();
        GetXaxis()->SetNoExponent();
        GetYaxis()->SetNdivisions(505);

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});

        GetFrame()->SetTitle("");
        if(i == 1)
            GetYaxis()->SetTitle("#LT p_{T}^{PUPI}/p_{T}^{CHS}#GT");
        if(i == 42)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} [GeV]");

        DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u));

        if(i == 42) {
            TLegend * leg = new TLegend(0.05, 0.45, 0.57, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            if(plotType == "time") {
                for(int k = 0; k < nRuns; ++k) {
                    TString t = k != 0 ? SF("Run %c",'A'+ k) : "All Runs";
                    leg->AddEntry(prof[k], t, "l");
                }
            }
            else if(plotType == "pileUp"){
                leg->AddEntry(prof[0], "any #rho", "l");
                leg->AddEntry(prof[1], "#rho < 15", "l");
                leg->AddEntry(prof[2], "15 < #rho < 22", "l");
                leg->AddEntry(prof[3], "22 < #rho < 30", "l");
                leg->AddEntry(prof[4], "30 < #rho", "l");
            }



            if(plotFit) {
                TF1 *funR = new TF1(SF("%d",rand()), "x");
                TF1 *funB = new TF1(SF("%d",rand()), "x");
                funB->SetLineColor(kBlue);
                leg->AddEntry(prof[0], "All Runs", "l");
                leg->AddEntry(funR, "LogLin fit", "l");
                leg->AddEntry(funB, "Extrap.", "l");
            }

            leg->Draw();
        }

        UpdateFrame();

    }
    
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
    return make_pair(hFit84, hFit1000);
}


void PLOTTER::MatchingFactorsYbased(vector<double> res, bool allPer)
{
    int nRuns = 1;
    if(allPer) nRuns = 8;


    //cout << "x01234 " << res[0] << " "<< res[1] << " "<< res[2] << " "<< res[3] << " "<< res[4] << " "<< endl;


    TH3D *hBalEtaPt;
    hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAll[0][0]->Clone(SF("%d",rand())));
    assert(hBalEtaPt);


    TCanvas *can = new TCanvas("can", "can", 800, 500);

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(0);

    //can->
    //dnPad->Divide(5, 4, 0.0001, 0.002);
    dnPad->cd();
    SetLeftRight(0.08, 0.03);
    DividePad(vector<double>(7,1.), vector<double>(6,1.));
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 7*6; ++i) {
        dnPad->cd(i);
        //gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
        //gPad->SetLogy();
        //if(i >= 19) continue;
        int iS = 19+1*i ;
        int iE = 19+1*(i+1) -1 ;


        double l = hBalEtaPt->GetYaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetYaxis()->GetBinUpEdge(iE);

        //cout << "xOrg " << hBalEtaPt->GetXaxis()->GetXmin() <<" "<< hBalEtaPt->GetXaxis()->GetXmax()<< endl;
        //cout << "yOrg " << hBalEtaPt->GetYaxis()->GetXmin() <<" "<< hBalEtaPt->GetYaxis()->GetXmax()<< endl;
        //cout << "zOrg " << hBalEtaPt->GetZaxis()->GetXmin() <<" "<< hBalEtaPt->GetZaxis()->GetXmax()<< endl;


        vector<TProfile*> prof(8);
        for(int k = 0; k < nRuns; ++k) {
            //iS = min(19,iS);
            //iE = min(19,iE);
            //cout << "RANGESstart " <<  hBalEtaPtAll[k][0]->GetXaxis()->GetXmin() << " "<< hBalEtaPtAll[k][0]->GetXaxis()->GetXmax() << endl;
            hBalEtaPtAll[0][k]->GetYaxis()->SetRange(iS,iE);
            hBalEtaPtAll[0][k]->GetXaxis()->SetRange(1,42);
            TH2D *hTemp = dynamic_cast<TH2D*>(hBalEtaPtAll[0][k]->Project3D(SF("%d_xz",rand()))); 
            assert(hTemp);
            prof[k] = hTemp->ProfileY();

            CleanProfile(prof[k], 100);

            //cout << "RANGES " << prof[k]->GetXaxis()->GetXmin() << " "<< prof[k]->GetXaxis()->GetXmax() << endl;

            prof[k]->SetLineColor(1+k);
            if(k==0)
                prof[k]->Draw("hist e");
            else
                prof[k]->Draw("hist e same");
        }
        prof[0]->Draw("hist e same");


        double pt = hBalEtaPt->GetYaxis()->GetBinCenter(iS);
        if(res.size() > 0) {
            //fun->SetParameter(0, pt);
            //fun->DrawCopy("same");
        }



        /*
        //Plot result of the fit
        TH1D *hfit1 = dynamic_cast<TH1D*>(hFit->ProjectionX(SF("%d",rand()), iS, iE));
        cout << "I am here " << __LINE__ << prof[0] << endl;
        hfit1->SetLineColor(kRed);
        cout << "I am here " << __LINE__ << prof[0] << endl;
        hfit1->Draw("hist same");
        */



        GetXaxis()->SetRangeUser(0, 5);
        GetYaxis()->SetRangeUser(0.961, 1.119);
        GetYaxis()->SetNdivisions(505);

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});

        GetFrame()->SetTitle("");
        if(i == 1)
            GetYaxis()->SetTitle("#LT p_{T}^{PUPI}/p_{T}^{CHS}#GT");
        if(i == 42)
            GetXaxis()->SetTitle("y");

        DrawLatexUp(-1, SF("%3.0f < p_{T} < %3.0f", l, u), 10);

        if(i == 42) {
            TLegend * leg = new TLegend(0.05, 0.45, 0.37, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            for(int k = 0; k < nRuns; ++k) {
                TString t = k != 0 ? SF("Run %c",'A'+ k) : "All Runs";
                leg->AddEntry(prof[k], t, "l");
            }
            //leg->AddEntry(prof[0], "any #rho", "l");
            //leg->AddEntry(prof[1], "any #rho", "l");
            //leg->AddEntry(prof[2], "any #rho", "l");
            //leg->AddEntry(prof[3], "any #rho", "l");


            //leg->AddEntry(prof[1], "#rho < 15", "l");
            //leg->AddEntry(prof[2], "15 < #rho < 22", "l");
            //leg->AddEntry(prof[3], "22 < #rho < 30", "l");
            //leg->AddEntry(prof[4], "30 < #rho", "l");
            leg->Draw();
        }



    }
    
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}











void PLOTTER::AsymmetryEtaDep()
{
    auto hBalEtaPt = hBalEtaPtAll[0][perID];
    auto hBalEtaPtUp = hBalEtaPtAllUp[perID];
    auto hBalEtaPtDn = hBalEtaPtAllDn[perID];

    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //can->
    dnPad->Divide(5, 4, 0.0001, 0.002);
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 20; ++i) {
        dnPad->cd(i);
        gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
        if(i >= 19) continue;
        int iS = 1*i ;
        int iE = 1*(i+1) -1 ;
        TH1D *hEta = hBalEtaPt->ProjectionZ(SF("ProjEta%d%d",i,rand()),  iS, iE, 0, -1 );
        TH1D *hEtaUp = hBalEtaPtUp->ProjectionZ(SF("ProjEtaUp%d%d",i,rand()),  iS, iE, 0, -1 );
        TH1D *hEtaDn = hBalEtaPtDn->ProjectionZ(SF("ProjEtaDn%d%d",i,rand()),  iS, iE, 0, -1 );

        double l = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);
        hEta->SetTitle(SF("%1.1f < |#eta| < %1.1f", l, u));
        hEta->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hEta->SetLineColor(kBlack);

        hEtaUp->SetLineColor(kBlue);
        hEtaDn->SetLineColor(kBlue);
        hEtaUp->SetLineStyle(2);
        hEtaDn->SetLineStyle(2);

        hEta->SetMaximum(1.2*hEta->GetMaximum());

        SetFontSizes(hEta, 0.07);
        hEta->Draw("hist e ");
        hEtaUp->Draw("hist e same");
        hEtaDn->Draw("hist e same");
        GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
    }
    
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}


void PLOTTER::AsymmetryPtDep()
{
    auto hBalEtaPt = hBalEtaPtAll[0][perID];
    auto hBalEtaPtUp = hBalEtaPtAllUp[perID];
    auto hBalEtaPtDn = hBalEtaPtAllDn[perID];

    TCanvas *can = new TCanvas("can", "can");
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    dnPad->Divide(5, 3, 0.0001, 0.002);

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 15; ++i) {
        dnPad->cd(i);
        gPad->SetMargin(0.14, 0.05, 0.18, 0.09);
        int iS = 3*i + 6;
        int iE = 3*(i+1) -1 + 6;
        TH1D *hPt = hBalEtaPt->ProjectionZ(SF("ProjPt%d%d",i,rand()),  0, -1, iS, iE );
        TH1D *hPtUp = hBalEtaPtUp->ProjectionZ(SF("ProjPtUp%d%d",i,rand()),  0, -1, iS, iE );
        TH1D *hPtDn = hBalEtaPtDn->ProjectionZ(SF("ProjPtDn%d%d",i,rand()),  0, -1, iS, iE );

        double l = hBalEtaPt->GetYaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetYaxis()->GetBinUpEdge(iE);
        hPt->SetTitle(SF("%3.0f < p_{T} < %3.0f", l, u));
        hPt->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hPt->SetMaximum(1.2*hPt->GetMaximum());
        SetFontSizes(hPt, 0.07);
        hPt->Draw("hist e");

        hPtUp->SetLineColor(kBlue);
        hPtDn->SetLineColor(kBlue);
        hPtUp->SetLineStyle(2);
        hPtDn->SetLineStyle(2);
        hPtUp->Draw("hist e same");
        hPtDn->Draw("hist e same");

    }
    //can->Print(outName +"(");
    can->Print(outName);
    can->Clear();
    delete can;
}





void PLOTTER::AsymmetryEtaPtDep(int eta1, int eta2)
{
    auto hBalEtaPt = hBalEtaPtAll[0][perID];
    auto hBalEtaPtUp = hBalEtaPtAllUp[perID];
    auto hBalEtaPtDn = hBalEtaPtAllDn[perID];

    vector<int> etaBins;
    for(int e = eta1; e <= eta2; ++e)
        etaBins.push_back(e);
    etaBins.push_back(eta2+1);
    int nEta = etaBins.size()-1;
    
    //k= {1, 4, 7, 10, 13, 16};
    const vector<int> ptBins = {12, 18, 24, 30, 36, 43, 46};
    int nPt = ptBins.size() - 1;

    //TProfile2D * hProf2 = hBalEtaPt->Project3DProfile();
    //Print grid

    TCanvas *can = new TCanvas("can", "can");
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //dnPad->Divide(5, 5, 0, 0);


    //dnPad->SetBottomMargin(   0.1);
    dnPad->cd();
    SetTopBottom(0.1, 0.1 + 0.8/5.*(5-nEta));

    DividePad(vector<double>(ptBins.size()-1,1.), vector<double>(nEta,1.));

    for(int i = 1; i <= nEta; ++i) 
    for(int j = 1; j <= nPt;  ++j) {
        dnPad->cd((i-1)*nPt + j);
        int iS = etaBins[i-1];
        int iE = etaBins[i+1-1]-1;
        int jS = ptBins[j-1];
        int jE = ptBins[j+1-1]-1;

        TH1D *hBoth   = hBalEtaPt->ProjectionZ(SF("ProjBoth%d%d%d",i,j,rand()),  iS, iE, jS, jE);
        TH1D *hBothUp = hBalEtaPtUp->ProjectionZ(SF("ProjBothUp%d%d%d",i,j,rand()),  iS, iE, jS, jE);
        TH1D *hBothDn = hBalEtaPtDn->ProjectionZ(SF("ProjBothDn%d%d%d",i,j,rand()),  iS, iE, jS, jE);


        hBoth->SetMaximum(1.2*hBoth->GetMaximum());
        //SetFontSizes(hBoth, 0.07);
        //hBoth->SetTitleSize(0.13);
        hBoth->Draw("hist e ");

        hBothUp->SetLineColor(kBlue);
        hBothDn->SetLineColor(kBlue);
        hBothUp->SetLineStyle(2);
        hBothDn->SetLineStyle(2);


        //hBothUp->Draw("hist e same");
        //hBothDn->Draw("hist e same");

        GetFrame()->SetTitle("");

        SetFTO({11}, {5.1}, {1.3, 3.0, 0.3, 1.1});
        GetYaxis()->SetLabelSize(0);
        GetXaxis()->SetNdivisions(505);
        GetYaxis()->SetNdivisions(505);

        if(j % nPt == 0)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");


        double lPt = hBalEtaPt->GetYaxis()->GetBinLowEdge(jS);
        double uPt = hBalEtaPt->GetYaxis()->GetBinUpEdge(jE);
        double lEt = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double uEt = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);
        //hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));

        if(i == 1)
            DrawLatexUp(dnPad->GetPad((i-1)*nPt + j),  1.1, SF("%3.0f < p_{T} < %3.0f", lPt, uPt));
        if(j == 1)
            DrawLatexLeft(dnPad->GetPad((i-1)*nPt + j),  2.1, SF("%1.1f < |#eta| < %1.1f", lEt, uEt), -1, ">" );

    }

    can->Print(outName);
    can->Clear();
    delete can;
}


void PLOTTER::AsymmetryEtaPtTimeDep(int eta1, int eta2)
{
    auto hBalEtaPt = hBalEtaPtAll[0][0];

    vector<int> etaBins;
    for(int e = eta1; e <= eta2; ++e)
        etaBins.push_back(e);
    etaBins.push_back(eta2+1);
    int nEta = etaBins.size()-1;
    
    //k= {1, 4, 7, 10, 13, 16};
    const vector<int> ptBins = {20, 24, 27, 30, 33, 36, 40, 43, 46, 49, 54, 58};
    int nPt = ptBins.size() - 1;

    //TProfile2D * hProf2 = hBalEtaPt->Project3DProfile();
    //Print grid

    TCanvas *can = new TCanvas("can", "can", 800, 500);
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //dnPad->Divide(5, 5, 0, 0);


    //dnPad->SetBottomMargin(   0.1);
    dnPad->cd();
    SetTopBottom(0.1, 0.1 + 0.8/5.*(5-nEta));
    SetLeftRight(0.05, 0.05);

    DividePad(vector<double>(ptBins.size()-1,1.), vector<double>(nEta,1.));

    for(int i = 1; i <= nEta; ++i) 
    for(int j = 1; j <= nPt;  ++j) {
        dnPad->cd((i-1)*nPt + j);
        int iS = etaBins[i-1];
        int iE = etaBins[i+1-1]-1;
        int jS = ptBins[j-1];
        int jE = ptBins[j+1-1]-1;

        TH1D *hBoth  = hBalEtaPtAll[0][0]->ProjectionZ(SF("ProjBoth%d%d%d",i,j,rand()),  iS, iE, jS, jE);
        hBoth->SetLineColor(kBlack);
        hBoth->Rebin(2);
        hBoth->Draw("hist e ");

        double Max = hBoth->GetMaximum();

        vector<TH1D*> hBothNow(8);
        hBothNow[0] = hBoth;
        for(int k = 1; k < 8; ++k) {
            hBothNow[k]  = hBalEtaPtAll[0][k]->ProjectionZ(SF("ProjBoth%d%d%d",i,j,rand()),  iS, iE, jS, jE);

            hBothNow[k]->Rebin(2);

            if(hBothNow[k]->Integral() > 0)
                hBothNow[k]->Scale(hBoth->Integral() / hBothNow[k]->Integral());
            hBothNow[k]->SetLineColor(k+1);
            hBothNow[k]->Draw("hist e same");
            Max = max(Max, hBothNow[k]->GetMaximum());
        }
        hBoth->Draw("hist e same");

        GetXaxis()->SetRangeUser(0.87, 1.13);
        GetFrame()->SetMaximum(1.2*Max);

        GetFrame()->SetTitle("");

        SetFTO({11}, {5.1}, {1.3, 3.0, 0.3, 1.1});
        GetYaxis()->SetLabelSize(0);
        GetXaxis()->SetNdivisions(505);
        GetYaxis()->SetNdivisions(505);

        if(j % nPt == 0)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");


        double lPt = hBalEtaPt->GetYaxis()->GetBinLowEdge(jS);
        double uPt = hBalEtaPt->GetYaxis()->GetBinUpEdge(jE);
        double lEt = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double uEt = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);
        //hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));

        if(i == 1)
            DrawLatexUp(dnPad->GetPad((i-1)*nPt + j),  1.1, SF("%3.0f < p_{T} < %3.0f", lPt, uPt), 8);
        if(j == 1)
            DrawLatexLeft(dnPad->GetPad((i-1)*nPt + j),  2.1, SF("%1.2f < |#eta| < %1.2f", lEt, uEt), 10, ">" );

        if(i == 1 && j == nPt) {
            TLegend * leg = new TLegend(0.05, 0.15, 0.32, 0.55);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            for(int k = 0; k < 8; ++k) {
                TString t = k != 0 ? SF("Run %c",'A'+ k) : "All Runs";
                leg->AddEntry(hBothNow[k], t, "l");
            }
            leg->Draw();

        }
    }

    can->Print(outName);
    can->Clear();
    delete can;
}

void PLOTTER::AsymmetryEtaPtPileUpDep(int eta1, int eta2)
{
    auto hBalEtaPt = hBalEtaPtAll[0][0];

    vector<int> etaBins;
    for(int e = eta1; e <= eta2; ++e)
        etaBins.push_back(e);
    etaBins.push_back(eta2+1);
    int nEta = etaBins.size()-1;
    
    //k= {1, 4, 7, 10, 13, 16};
    const vector<int> ptBins = {12, 18, 24, 30, 36, 43, 46};
    int nPt = ptBins.size() - 1;

    //TProfile2D * hProf2 = hBalEtaPt->Project3DProfile();
    //Print grid

    TCanvas *can = new TCanvas("can", "can");
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //dnPad->Divide(5, 5, 0, 0);


    //dnPad->SetBottomMargin(   0.1);
    dnPad->cd();
    SetTopBottom(0.1, 0.1 + 0.8/5.*(5-nEta));

    DividePad(vector<double>(ptBins.size()-1,1.), vector<double>(nEta,1.));

    for(int i = 1; i <= nEta; ++i) 
    for(int j = 1; j <= nPt;  ++j) {
        dnPad->cd((i-1)*nPt + j);
        int iS = etaBins[i-1];
        int iE = etaBins[i+1-1]-1;
        int jS = ptBins[j-1];
        int jE = ptBins[j+1-1]-1;

        TH1D *hBoth  = hBalEtaPtAll[0][0]->ProjectionZ(SF("ProjBoth%d%d%d",i,j,rand()),  iS, iE, jS, jE);
        hBoth->SetLineColor(kBlack);

        hBoth->Rebin(4);

        hBoth->Draw("hist e ");

        double Max = hBoth->GetMaximum();

        vector<TH1D*> hBothNow(5);
        hBothNow[0] = hBoth;
        for(int k = 1; k < 5; ++k) {
            hBothNow[k]  = hBalEtaPtAll[k][0]->ProjectionZ(SF("ProjBoth%d%d%d",i,j,rand()),  iS, iE, jS, jE);

            hBothNow[k]->Rebin(4);

            if(hBothNow[k]->Integral() > 0)
                hBothNow[k]->Scale(hBoth->Integral() / hBothNow[k]->Integral());
            hBothNow[k]->SetLineColor(k+1);
            hBothNow[k]->Draw("hist e same");
            Max = max(Max, hBothNow[k]->GetMaximum());
        }
        hBoth->Draw("hist e same");

        GetFrame()->SetMaximum(1.2*Max);

        GetFrame()->SetTitle("");

        SetFTO({11}, {5.1}, {1.3, 3.0, 0.3, 1.1});
        GetYaxis()->SetLabelSize(0);
        GetXaxis()->SetNdivisions(505);
        GetYaxis()->SetNdivisions(505);

        if(j % nPt == 0)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");


        double lPt = hBalEtaPt->GetYaxis()->GetBinLowEdge(jS);
        double uPt = hBalEtaPt->GetYaxis()->GetBinUpEdge(jE);
        double lEt = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double uEt = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);
        //hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));

        if(i == 1)
            DrawLatexUp(dnPad->GetPad((i-1)*nPt + j),  1.1, SF("%3.0f < p_{T} < %3.0f", lPt, uPt));
        if(j == 1)
            DrawLatexLeft(dnPad->GetPad((i-1)*nPt + j),  2.1, SF("%1.1f < |#eta| < %1.1f", lEt, uEt), -1, ">" );

        if(i == 1 && j == nPt) {
            TLegend * leg = new TLegend(0.05, 0.15, 0.32, 0.55);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            //for(int k = 0; k < 5; ++k) {
                //TString t = k != 0 ? SF("Run %c",'A'+ k) : "All Runs";
                //leg->AddEntry(hBothNow[k], t, "l");
            //}

            leg->AddEntry(hBothNow[0], "any #rho", "l");
            leg->AddEntry(hBothNow[1], "#rho < 15", "l");
            leg->AddEntry(hBothNow[2], "15 < #rho < 22", "l");
            leg->AddEntry(hBothNow[3], "22 < #rho < 30", "l");
            leg->AddEntry(hBothNow[4], "30 < #rho", "l");
            leg->Draw();


            leg->Draw();

        }
    }

    can->Print(outName);
    can->Clear();
    delete can;
}



Double_t fun2D(Double_t *x, Double_t *par) {
    Double_t lx = log(x[0]);

    Double_t p0 = par[0];

    //Double_t p1 = par[1];
    //Double_t p2 = par[2];

    Double_t a2 = par[1];
    Double_t a3 = par[2];
    Double_t b2 = par[3];
    Double_t b3 = par[4];


    //Double_t result = p0 + p1*lx + p2*x[1];
    Double_t result = p0 + (a2*lx + a3*lx*lx) + (b2*x[1] + b3*x[1]*x[1]);

    return result;
}


double polFun(double pt, double y, const double *x)
{
    double lpt = log10(pt);

    double fun  = x[0] + x[1]*lpt + x[2]*y;
           fun += x[3]*lpt*lpt + x[4]*y*y;
    return fun;
}


class  MyChi2 {
    public:
        TProfile2D *prof;
        double yMax=3.5, pMin=80, pMax=900;
        int npar = 5;
        // use constructor to customize your function object
        double operator() (const double *x) {
            // function implementation using class data members

            //return (pow(x[0]-3,2) + pow(x[1]-2,2));

            double sum2 = 0;

            for(int i = 0; i <= prof->GetNbinsX(); ++i)
            for(int j = 0; j <= prof->GetNbinsY(); ++j) {
                double y  = prof->GetYaxis()->GetBinCenter(j);
                double pt = prof->GetXaxis()->GetBinCenter(i);

                if(y < 0 ||  y > yMax || pt < pMin || pt > pMax) continue;
                double val = prof->GetBinContent(i,j);
                double err = prof->GetBinError(i,j);

                double fun = polFun(pt, y, x);

                if(err > 0)
                    sum2 += pow((val - fun)/err, 2);
            }
            //cout << sum2 << endl;
            return sum2;

        }
};


vector<double> NumericalMinimization(TProfile2D *prof)
{
    // Choose method upon creation between:
    // kConjugateFR, kConjugatePR, kVectorBFGS,
    // kVectorBFGS2, kSteepestDescent
    ROOT::Math::GSLMinimizer min( ROOT::Math::kVectorBFGS );

    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(300000);
    min.SetTolerance(0.001);

    MyChi2 mychi2;
    mychi2.prof = prof;


    ROOT::Math::Functor f(mychi2, mychi2.npar); 
    vector<double> step = {0.01, 0.01, 0.01, 0.01, 0.01};
    vector<double> variable = {1.10669,-0.0340574, 0.00321177, 0, 0};

    min.SetFunction(f);

    // Set the free variables to be minimized!

    for(int i = 0; i < mychi2.npar; ++i)
        min.SetVariable(i, SF("x%d",i).Data(), variable[i], step[i]);

    //min.SetVariable(0,"x",variable[0], step[0]);
    //min.SetVariable(1,"y",variable[1], step[1]);

    cout << "Start min" << endl;
    min.Minimize(); 
    cout << "End min" << endl;

    const double *xs = min.X();
    vector<double> res(xs, xs + mychi2.npar);

    //cout << "Minimum: f(" << xs[0] << "," << xs[1] <<" "<< xs[2] <<  "): " 
        //<< mychi2(xs) << endl;

    for(int i = 0; i <= prof->GetNbinsX(); ++i)
    for(int j = 0; j <= prof->GetNbinsY(); ++j) {
        double y  = prof->GetYaxis()->GetBinCenter(j);
        double pt = prof->GetXaxis()->GetBinCenter(i);

        if(y < 0 ||  y > mychi2.yMax || pt < mychi2.pMin || pt > mychi2.pMax) continue;
        double val = prof->GetBinContent(i,j);
        double err = prof->GetBinError(i,j);

        double fun = polFun(pt, y, xs);

        //cout << "pt,y " << pt << " " << y << " : " << fun << " "<< val << endl;
    }
    cout << "x01234 " << xs[0] << " "<< xs[1] << " "<< xs[2] << " "<< xs[3] << " "<< xs[4] << " "<< endl;

    return res;

}








vector<double> PLOTTER::MeanAsym(int shift, TString style)
{

    TH3D *hBalEtaPt;
    if(shift == 0)
        hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAll[0][perID]->Clone(SF("%d",rand())));
    else if(shift == 1)
        hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAllUp[perID]->Clone(SF("%d",rand())));
    else
        hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAllDn[perID]->Clone(SF("%d",rand())));

    assert(hBalEtaPt);

    //hBalEtaPt->RebinY(4);
    gStyle->SetPaintTextFormat("2.2f");
    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);


    dnPad->SetMargin(0.1, 0.1, 0.1, 0.1);
    dnPad->SetLogx();

    TProfile2D * hProf = hBalEtaPt->Project3DProfile();
    hProf->GetXaxis()->SetRangeUser(84, 2000);

    hProf->BuildOptions(0.8, 1.2, "s");
    CleanProfile(hProf, 100);

    hProf->SetTitle("");
    /*
    if(shift == 0)
        hProf->SetTitle("Mean value of p_{T}^{PUPPI} / p_{T}^{CHS}");
    else if(shift == 1)
        hProf->SetTitle("Mean value of p_{T}^{PUPPI} / p_{T}^{CHS} (JECup)");
    else
        hProf->SetTitle("Mean value of p_{T}^{PUPPI} / p_{T}^{CHS} (JECdn)");
    */


    hProf->GetYaxis()->SetTitle("#eta");
    hProf->GetXaxis()->SetTitle("p_{T}");
    hProf->GetXaxis()->SetMoreLogLabels();

    hProf->SetMaximum(1.15);
    hProf->SetMinimum(0.85);


    //Do fit
    //cout << "Starting fitting" << endl;
    //vector<double> res = NumericalMinimization(hProf);

    /*

    f2 = new TF2("fun2D",fun2D, 84,150, 0, 0.4, 5);
    //f2->SetParameters(1, 0, 0);
    f2->SetParameters(1, 0, 0, 0, 0);


    hFit = (TH2D*) hProf->ProjectionXY("tada"); 
    hFit->GetXaxis()->SetRangeUser(84, 2000);

    hProf->Fit("fun2D", "eN");


    for(int i = 1; i < hFit->GetNbinsX(); ++i)
    for(int j = 1; j < hFit->GetNbinsY(); ++j) {
        double x = hFit->GetXaxis()->GetBinCenter(i);
        double y = hFit->GetYaxis()->GetBinCenter(j);

        double val = f2->Eval(x,y);
        cout << "x y val " << x <<" "<< y <<" " << val << endl;
        hFit->SetBinContent(i,j, val);
        hFit->SetBinError(i,j, 0);
    }
    */


    if(style == "colz")
        hProf->Draw("colz");
    else
        hProf->Draw("text");



    can->Print(outName);
    can->Clear();
    delete hProf;
    delete can;


    /*
    can = new TCanvas("can", "can");

    tie(upPad, dnPad) = TitleSpace(perID);


    dnPad->SetMargin(0.1, 0.1, 0.1, 0.1);
    dnPad->SetLogx();

    hFit->Draw("text");

    can->Print(outName);
    */

    return {0,0,0};
}







void plotter()
{

    TH1::SetDefaultSumw2();


    PLOTTER plot;
    //TString myOut = "plots/Spring16_25nsV6.pdf";
    //plot.Init("histos/Spring16_25nsV6__Spring16_25nsV6.root", myOut+"(");

    //TString myOut = "plots/Summer16_07Aug2017V5newnew.pdf";
    //plot.Init("histos/Summer16_07Aug2017V5__Summer16_07Aug2017V5newnew.root", myOut+"(");

    //TString myOut = "plots/new10.pdf";
    //plot.Init("histos/V10V10newBinning.root", myOut+"(");
    //TString myOut = "plots/new7.pdf";
    //plot.Init("histos/V7V7newBinning.root", myOut+"(");

    //TString nTag = "profilesCentralBB";
    TString nTag = "newAsym4ak8";


    TString myOut = "plots/"+nTag+".pdf";
    //plot.Init("histos/V11V11newBinning.root", myOut+"(");
    plot.Init("/nfs/dust/cms/user/zlebcr/JEC/ntuplesTemp/histos/histo4AK8CHS.root", myOut+"(");


    //plot.Init("histos/Summer16_07Aug2017V5__Summer16_07Aug2017V5noResNew.root", myOut+"(");

    //TString myOut = "plots/Summer16_07Aug2017V5noRes.pdf";
    //plot.Init("histos/Summer16_07Aug2017V5__Summer16_07Aug2017V5noRes.root", myOut+"(");


    //TString myOut = "plots/Mixed.pdf";
    //plot.Init("histos/Summer16_07Aug2017V5__Spring16_23Sep2016V2new.root", myOut+"(");
    //TString myOut = "plots/MixedNoRes.pdf";
    //plot.Init("histos/Summer16_07Aug2017V5__Spring16_23Sep2016V2noRes.root", myOut+"(");



    //plot.Init("histoDir/Aug/jetsAll.root", myOut+"(");
    /*
    plot.perID = 1;

    plot.JEC();
    plot.outName = myOut;
    plot.PtEtaDep();
    plot.outName = myOut+ ")";
    plot.MeanAsym(0);
    return;
    */
    gStyle->SetOptStat(0);
    //plot.PlotProfilesOverlap();
    //return;

    plot.perID = 0;
    plot.Unmatched();
    plot.outName = myOut;
    plot.Unmatched("chs4");
    plot.PlotProfiles(0, -1);
    plot.PlotProfiles(1, -1);
    plot.PlotProfiles(2, -1);
    plot.PlotProfiles(3, -1);

    plot.PlotProfiles(7, -1);
    plot.PlottMatchingCorr(+1, false);
    plot.PlottMatchingCorr(-1, false);
    TH1D *h1P, *h2P, *h1N, *h2N;
    tie(h1P, h2P) = plot.PlottMatchingCorr(+1, true);

    tie(h1N, h2N) = plot.PlottMatchingCorr(-1, true);
    //plot.PlotFit(h1P, h2P);
    plot.PlotFitFour(h1P, h2P, h1N, h2N);

    plot.PlotProfiles(5, -1);
    plot.PlotProfiles(6, -1);
    plot.PlotProfiles(8, -1);
    plot.PlotProfiles(9, -1);
    //plot.PlotProfiles(3, -1, true);
    //plot.PlotProfiles(7, -1, true);
    plot.PlotVar1overVar2();

    plot.CompareSigmas();

    plot.PlotProfiles(10, -1);

    plot.outName = myOut + ")";
    plot.PlotProfiles(4, -1);
    return;
    /*
    for(int i = 0; i < 8; ++i)
        plot.AsymmetryEtaPtTimeDep(5*i+1, 5*i+5);
    plot.AsymmetryEtaPtTimeDep(5*8+1, 5*8+2);
    */

    //plot.AsymmetryEtaPtTimeDep(6, 10);
    //plot.AsymmetryEtaPtTimeDep(11, 15);
    //plot.AsymmetryEtaPtTimeDep(16, 18);
    /*
    plot.MatchingFactorsTimeDep(0);
    //plot.MatchingFactorsTimeDep(1);
    //plot.MatchingFactorsTimeDep(2);
    //plot.MatchingFactorsTimeDep(3);


    plot.AsymmetryEtaPtPileUpDep(1, 5);
    plot.AsymmetryEtaPtPileUpDep(6, 10);
    plot.AsymmetryEtaPtPileUpDep(11, 15);
    plot.AsymmetryEtaPtPileUpDep(16, 18);
    */

    auto res = plot.MeanAsym(0, "colz");


    plot.MatchingFactorsPileUpDep(false,  "time");
    plot.MatchingFactorsPileUpDep(false,  "pileUp");

    TH1D *h1, *h2;
    tie(h1, h2) = plot.MatchingFactorsPileUpDep(true, "");
    plot.PlotFit(h1, h2);

    plot.MatchingFactorsYbased({}, true);
    plot.outName = myOut + ")";
    plot.MatchingFactorsYbased(res, false);

    //plot.MatchingFactorsTimeDep(4);
    return;



    int iMax = 7;
    for(int i = 0; i <= iMax; ++i) {
        plot.perID = i;
        plot.AsymmetryEtaDep();
        plot.outName = myOut;
        plot.AsymmetryPtDep();
        plot.AsymmetryEtaPtDep(1, 5);
        plot.AsymmetryEtaPtDep(6, 10);
        plot.AsymmetryEtaPtDep(11, 15);
        plot.AsymmetryEtaPtDep(16, 18);

        plot.PtEtaDep();
        plot.Unmatched();

        plot.MeanAsym(0);
        plot.MeanAsym(0, "colz");
        plot.MeanAsym(1);
        if(i == iMax)
            plot.outName = myOut + ")";
        plot.MeanAsym(2);
    }

    /*

    TH2D *hSigma = hProf->ProjectionXY("sigma");  

    for(int i = 1; i <= hProf->GetNbinsX(); ++i)
    for(int j = 1; j <= hProf->GetNbinsY(); ++j) {
        TH1D *hTemp = hBalEtaPt->ProjectionZ(SF("test%d%d",i,j), i, i, j, j);
        double s = 0;
        for(int k = 1; k <= hTemp->GetNbinsX(); ++k)
            s += pow(hTemp->GetBinCenter(k) - hTemp->GetMean(), 2) * hTemp->GetBinContent(k);
        s /= hTemp->Integral();

        cout << "Holka " << hTemp->GetStdDev()  << endl;
        hSigma->SetBinContent(i, j, hTemp->GetStdDev() );
        //hSigma->SetBinContent(i, j, sqrt(s) );
        delete hTemp;
    }

    hSigma->SetTitle("STD of p_{T}^{PUPPI} / p_{T}^{CHS}");
    hSigma->GetYaxis()->SetTitle("#eta");
    hSigma->GetXaxis()->SetTitle("p_{T}");

    hSigma->Draw("text");

    hJetPt->Draw();
    can->Print("res.pdf)");

    */
}
