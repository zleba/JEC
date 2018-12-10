#include "../bin/JEC.h"
#include <map>
#include <cmath>
#include "TString.h"


const vector<double> etaBins2 =  {0., 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.538, 4.716, 4.889, 5.191};


const vector<double> Ptbinning = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,  2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717,    7000};


//const vector<TString> pers  = {"MC", "BCD", "EF", "GH"};
const vector<TString> pers1  = {"BCD", "EF", "GH", "MC"};
const vector<TString> pers2  = {"B","C", "D", "E", "F", "MC"};
const vector<TString> types = {"L2L3res", "L2rel"};

const int coneSize = 4;
const int version1 = 18;
const int version2 = 28;

const string jecTag1 = "Summer16_07Aug2017";
//const string jecTag2 = "Summer16_07Aug2017";
//const string jecTag1 = "Fall17_17Nov2017";
const string jecTag2 = "Fall17_17Nov2017";


#include "/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/libPlottingHelper.so)

using namespace PlottingHelper;


using namespace std;


double GetPmax (double eta)
{
    //E = pT * cosh(eta);
    const double Ebeam = 6500;
    return Ebeam / cosh(eta);
}
/*
void CleanHistos(vector<TH1D*> hArr, eta) 
{
    double pTmax = GetPmax(eta);
    for(auto &h : hArr) {
        for(int i = 1; i <= h->GetNbinsX(); ++i) {
            double pt = h->GetXaxis()->GetBinCenter(i);
            if(pt >= pTmax)
                h->SetBinContent(i, 0);
        }
    }
}
*/

class JECPlotter {
    public:


        struct JECset {
            TH2D *h1p, *h2p;
            TH2D *h1n, *h2n;
        };
        map<TString, map<TString, JECset>>  hJEC1, hJEC2;

        TString outFile;

        JECs *InitJEC(TString jecTag, int version, TString period) {

            //bool isMC = period == "M" ? true : false;
            //vector<string> dumy;

            TString type = (period == "MC") ? "MC" : "DATA";
            if(type == "MC") period = "";

            //jecMapCHS4[per]   = new JECs("Summer16_07Aug2017BCD_V14_DATA", "AK4PFchs");
            //jecMapCHS8[per]   = new JECs("Summer16_07Aug2017BCD_V14_DATA", "AK8PFchs");
            //jecMapPUPPI4[per] = new JECs("Summer16_07Aug2017BCD_V14_DATA", "AK4PFPuppi");

            JECs *jetCorrs = new JECs(Form("%s%s_V%d_%s",jecTag.Data(), period.Data(), version, type.Data()), Form("AK%dPFchs",coneSize));
            //jetCorrs2 = new JECs(Form("%s%s_V%d_%s",jecTag2.c_str(), period.Data(), version2, type.Data()), Form("AK%dPFchs",coneSize));

            //jetCorrsCHS.Init(isMC, jecTagCHS, period, versionCHS, SF("AK%dPFchs", coneSize).Data(), "", dumy);
            //jetCorrsPUPPI.Init(isMC, jecTagCHS, period, versionCHS, SF("AK%dPFPuppi", coneSize).Data(), "", dumy);
            return jetCorrs;
        }
        map<TString, map<TString, JECset>>  InitHistos(vector<TString> types, const vector<TString> &pers) {
            map<TString, map<TString, JECset>> hJEC;
            
            for(auto t : types) 
            for(auto p : pers)  {
                TString title = t + " (" + p + ")";
                hJEC[t][p].h1p = new TH2D(p+t+SF("%dp",rand()), title, etaBins2.size()-1, etaBins2.data(),
                                                                      Ptbinning.size()-1, Ptbinning.data());
                hJEC[t][p].h1n = new TH2D(p+t+SF("%dn",rand()), title, etaBins2.size()-1, etaBins2.data(),
                                                                      Ptbinning.size()-1, Ptbinning.data());
                assert(hJEC[t][p].h1p);
                assert(hJEC[t][p].h1n);
            }
            return hJEC;
        }


        void FillCorrs() {

            //funArr["L3abs"] = &JECs::ApplyJECL3Absolute ;
            hJEC1 = InitHistos(types, pers1);
            hJEC2 = InitHistos(types, pers2);
            
            FillCorr(hJEC1, types, pers1, jecTag1, version1);
            FillCorr(hJEC2, types, pers2, jecTag2, version2);
        }

        void FillCorr(map<TString, map<TString, JECset>> &hJEC,  vector<TString> types, vector<TString> pers, TString jecTag, int version)
        {
            typedef double (JECs::*FunPtr)(double &pt, double eta);
            map<TString, FunPtr > funArr;
            funArr["L2L3res"] = &JECs::ApplyJECL2L3Residual;
            funArr["L2rel"] = &JECs::ApplyJECL2Relative ;
            for(auto t : types) 
            for(auto p : pers) {
                //char per = p[0];

                //JECs *jetCorrs1, *jetCorrs2;
                JECs *jetCorrs =  InitJEC(jecTag, version, p);
                //InitJEC(jetCorrs2, p);

                //Fill Corrections
                //cout << "CorrHist " << hJEC[t][p].h1p << endl;
                for(int i = 1; i <= hJEC[t][p].h1p->GetNbinsX(); ++i)
                for(int j = 1; j <= hJEC[t][p].h1p->GetNbinsY(); ++j) {
                     double eta = hJEC[t][p].h1p->GetXaxis()->GetBinCenter(i);
                     double pt  = hJEC[t][p].h1p->GetYaxis()->GetBinCenter(j);

                     double pTmax =  GetPmax(eta);
                     
                     auto FillBin = [&](TH2D *h, double corr) {
                         if(pt < pTmax) {
                             h->SetBinContent(i, j, corr);
                             h->SetBinError(i, j, 0);
                         }
                     };

                     if(t != "L2L3res" || p != "MC")
                     {
                         double pt1,pt2;
                         pt1=pt2= pt;
                         FillBin(hJEC[t][p].h1p, ((*jetCorrs).*(funArr[t]))(pt1, eta));
                         FillBin(hJEC[t][p].h1n, ((*jetCorrs).*(funArr[t]))(pt2,-eta));

                         //FillBin(hJEC[t][p].h2p, ((*jetCorrs2).*(funArr[t]))(pt3, eta));
                         //FillBin(hJEC[t][p].h2n, ((*jetCorrs2).*(funArr[t]))(pt4,-eta));
                     }

                }
            }
        }




        void PlotAll() {

            const TString OutFile = SF("plots/jecPlotAK%d_v%d_%d.pdf",coneSize,version1,version2);

            auto plotOne = [&](TString t, TString p, double Min, double Max ) {
                //cout << "Radek " << hJEC[t][p].h1p << endl;
                PlottSingle(hJEC1[t][p].h1p,  hJEC1[t][p].h1n,
                            hJEC2[t][p].h1p, hJEC2[t][p].h1n, Min, Max);
            };

            double Min, Max;

            Min = 0.87, Max = 1.19;
            outFile = OutFile + "";
            //plotOne("L2L3res", "MC", Min, Max);
            //2016

            TString t = "L2L3res";
            PlottSingle(hJEC1[t]["GH"].h1p,  hJEC1[t]["GH"].h1n,
                            hJEC2[t]["F"].h1p, hJEC2[t]["F"].h1n, Min, Max);


            /*
            for(int p = 0; p < pers.size()-1; ++p) {
                plotOne("L2L3res", pers[p], Min, Max);
                if(p == 0) outFile = OutFile + "";
            }



            Min = 0.97, Max = 1.37;
            for(int p = 0; p < pers.size(); ++p) {
                if(p == pers.size() -1) outFile = OutFile + ")";
                plotOne("L2rel", pers[p], Min, Max);
            }

            */

            //plotOne("L2rel", "EF", Min, Max);
            //plotOne("L2rel", "GH", Min, Max);
            //plotOne("L2rel", "MC", Min, Max);

            /*
            Min = 0.96, Max = 1.04;
            plotOne("L3abs", "BCD", Min, Max);
            plotOne("L3abs", "EF", Min, Max);
            plotOne("L3abs", "GH", Min, Max);
            outFile = OutFile + ")";
            plotOne("L3abs", "MC", Min, Max);
            */


        }

        void PlottSingle(TH2D *h1p, TH2D *h1n, TH2D *h2p, TH2D *h2n, double Min=0.89, double Max=1.2) {
            TCanvas *can = new TCanvas("can", "can", 800, 500);

            SetLeftRight(0.08, 0.03);
            DividePad(vector<double>(7,1.), vector<double>(6,1.));

            gStyle->SetOptStat(0);

            for(int i = 1; i <= 7*6; ++i) {
                can->cd(i);
                //gPad->SetMargin(0.14, 0.05, 0.18, 0.08);
                gPad->SetLogx();
                //gPad->SetLogy();
                //if(i >= 19) continue;
                int iS = 1*i ;
                int iE = 1*(i+1) -1 ;


                double l = h1p->GetXaxis()->GetBinLowEdge(iS);
                double u = h1p->GetXaxis()->GetBinUpEdge(iE);
                double c = h1p->GetXaxis()->GetBinCenter(iE);

                TH1D *h1pProj   =  dynamic_cast<TH1D*>( h1p->ProjectionY(SF("%d",rand()), iS, iE) );
                TH1D *h2pProj =  dynamic_cast<TH1D*>( h2p->ProjectionY(SF("%d",rand()), iS, iE) );
                TH1D *h1nProj   =  dynamic_cast<TH1D*>( h1n->ProjectionY(SF("%d",rand()), iS, iE) );
                TH1D *h2nProj =  dynamic_cast<TH1D*>( h2n->ProjectionY(SF("%d",rand()), iS, iE) );

                h1pProj->Draw("hist e ][");
                h1pProj->SetLineColor(1);

                h2pProj->Draw("hist e same ][");
                h2pProj->SetLineColor(2);

                h1nProj->Draw("hist e same ][");
                h1nProj->SetLineColor(1);
                h1nProj->SetLineWidth(2);
                h1nProj->SetLineStyle(2);

                h2nProj->Draw("hist e same ][");
                h2nProj->SetLineColor(2);
                h2nProj->SetLineWidth(2);
                h2nProj->SetLineStyle(2);


                GetXaxis()->SetRangeUser(84, 2000);
                GetYaxis()->SetRangeUser(Min, Max);
                GetXaxis()->SetMoreLogLabels();
                GetXaxis()->SetNoExponent();
                GetYaxis()->SetNdivisions(505);


                //double pMax = GetPmax(c);
                //TLine *line = new TLine;
                //line->SetLineStyle(2);
                //line->DrawLine(pMax, Min, pMax, 1);


                SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 4.4});

                GetFrame()->SetTitle("");
                if(i == 1)
                    GetYaxis()->SetTitle("JEC");
                if(i == 42)
                    GetXaxis()->SetTitle("p_{T} [GeV]");

                DrawLatexUp(-1, SF("%1.2f < |#eta| < %1.2f", l, u));

                if(i == 42) {
                    TLegend * leg = new TLegend(0.15, 0.45, 0.87, 0.84);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->AddEntry(h1pProj, "2016JEC pos", "l");
                    leg->AddEntry(h1nProj, "2016JEC neg", "l");
                    leg->AddEntry(h2pProj, "2017JEC pos", "l");
                    leg->AddEntry(h2nProj, "2017JEC neg", "l");

                    leg->Draw();
                }

        }
        DrawLatexUp(can->GetPad(1), can->GetPad(7), 1.5, h1p->GetTitle() );
    
        can->Print(outFile);
        can->Clear();
        delete can;
    }

};

void JECplotter()
{

    JECPlotter jecPlotter;

    //jecPlotter.InitHistos();
    jecPlotter.FillCorrs();
    jecPlotter.PlotAll();




}
