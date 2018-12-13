#include "../bin/JEC.h"
#include <map>
#include <cmath>
#include "TString.h"


const vector<double> etaBins2 =  {0., 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.538, 4.716, 4.889, 5.191};


const vector<double> Ptbinning = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,  2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717,    7000};


//const vector<TString> pers  = {"MC", "BCD", "EF", "GH"};
const vector<TString> pers  = {"BCD", "EF", "GH", "MC"};
const vector<TString> types = {"L2L3res", "L2rel", "L3abs"};

const int coneSize = 4;
const int versionCHS   = 16;
const int versionPUPPI = 16;


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
            TH2D *hCHSp, *hPUPPIp;
            TH2D *hCHSn, *hPUPPIn;
        };

        map<TString, map<TString, JECset>>  hJEC;

        TString outFile;

        void InitJEC(JECs *&jetCorrsCHS, JECs *&jetCorrsPUPPI, TString  period) {
            const string jecTagCHS = "Summer16_07Aug2017";
            //bool isMC = period == "M" ? true : false;
            //vector<string> dumy;

            TString type = (period == "MC") ? "MC" : "DATA";
            if(type == "MC") period = "";
            if(type == "MC") {
                type = "DATA";
                period = "BCD";
            }

            //jecMapCHS4[per]   = new JECs("Summer16_07Aug2017BCD_V14_DATA", "AK4PFchs");
            //jecMapCHS8[per]   = new JECs("Summer16_07Aug2017BCD_V14_DATA", "AK8PFchs");
            //jecMapPUPPI4[per] = new JECs("Summer16_07Aug2017BCD_V14_DATA", "AK4PFPuppi");

            jetCorrsCHS   = new JECs(Form("Summer16_07Aug2017%s_V%d_%s",period.Data(),versionCHS,  type.Data()), Form("AK%dPFchs",coneSize));
            jetCorrsPUPPI = new JECs(Form("Summer16_07Aug2017%s_V%d_%s",period.Data(),versionPUPPI,type.Data()), Form("AK%dPFPuppi",coneSize));

            //jetCorrsCHS.Init(isMC, jecTagCHS, period, versionCHS, SF("AK%dPFchs", coneSize).Data(), "", dumy);
            //jetCorrsPUPPI.Init(isMC, jecTagCHS, period, versionCHS, SF("AK%dPFPuppi", coneSize).Data(), "", dumy);
        }
        void InitHistos() {

            for(auto t : types) 
            for(auto p : pers)  {

                TString title = t + " (" + p + ")";
                hJEC[t][p].hCHSp = new TH2D(p+t+"CHSp", title, etaBins2.size()-1, etaBins2.data(),
                                                                      Ptbinning.size()-1, Ptbinning.data());
                hJEC[t][p].hCHSn = new TH2D(p+t+"CHSn", title, etaBins2.size()-1, etaBins2.data(),
                                                                      Ptbinning.size()-1, Ptbinning.data());

                hJEC[t][p].hPUPPIp = new TH2D(p+t+"PUPPIp", title, etaBins2.size()-1, etaBins2.data(),
                                                                      Ptbinning.size()-1, Ptbinning.data());
                hJEC[t][p].hPUPPIn = new TH2D(p+t+"PUPPIn", title, etaBins2.size()-1, etaBins2.data(),
                                                                      Ptbinning.size()-1, Ptbinning.data());

            }
            //hCHSp["BCD"]   = new TH2D("hCHSp_BCD", "hCHSp_BCD", etaBins2.size()-1, etaBins2.data(), Ptbinning.size()-1, Ptbinning.data());
            //hPUPPIp["BCD"] = new TH2D("hPUPPIp_BCD", "hPUPPIp_BCD", etaBins2.size()-1, etaBins2.data(), Ptbinning.size()-1, Ptbinning.data());
        }

        void FillCorr() {

             typedef double (JECs::*FunPtr)(double &pt, double eta);

             map<TString, FunPtr > funArr;
             funArr["L2L3res"] = &JECs::ApplyJECL2L3Residual;
             funArr["L2rel"] = &JECs::ApplyJECL2Relative ;
             funArr["L3abs"] = &JECs::ApplyJECL3Absolute ;

            for(auto t : types) 
            for(auto p : pers) {
                char per = p[0];

                JECs *jetCorrsCHS, *jetCorrsPUPPI;
                InitJEC(jetCorrsCHS, jetCorrsPUPPI, p);

                //Fill Corrections
                for(int i = 1; i <= hJEC[t][p].hCHSp->GetNbinsX(); ++i)
                for(int j = 1; j <= hJEC[t][p].hCHSp->GetNbinsY(); ++j) {
                     double eta = hJEC[t][p].hCHSp->GetXaxis()->GetBinCenter(i);
                     double pt  = hJEC[t][p].hCHSp->GetYaxis()->GetBinCenter(j);
                     double cCHS, cPUPPI;

                     double pTmax =  GetPmax(eta);
                     
                     auto FillBin = [&](TH2D *h, double corr) {
                         if(pt < pTmax) {
                             h->SetBinContent(i, j, corr);
                             h->SetBinError(i, j, 0);
                         }
                     };

                     if(t != "L2L3res" || p != "MC")
                     {
                         double pt1,pt2,pt3,pt4;
                         pt1=pt2=pt3=pt4 = pt;
                         FillBin(hJEC[t][p].hCHSp, ((*jetCorrsCHS).*(funArr[t]))(pt1, eta));
                         FillBin(hJEC[t][p].hCHSn, ((*jetCorrsCHS).*(funArr[t]))(pt2,-eta));

                         FillBin(hJEC[t][p].hPUPPIp, ((*jetCorrsPUPPI).*(funArr[t]))(pt3, eta));
                         FillBin(hJEC[t][p].hPUPPIn, ((*jetCorrsPUPPI).*(funArr[t]))(pt4,-eta));
                     }

                }
            }
        }

        void PlotAll() {

            const TString OutFile = SF("plots/jecPlotAK%d_v%d.pdf",coneSize,versionCHS);

            auto plotOne = [&](TString t, TString p, double Min, double Max ) {
                PlottSingle(hJEC[t][p].hCHSp,  hJEC[t][p].hCHSn,
                            hJEC[t][p].hPUPPIp, hJEC[t][p].hPUPPIn, Min, Max);
            };

            double Min, Max;

            Min = 0.87, Max = 1.19;
            outFile = OutFile + "(";
            //plotOne("L2L3res", "MC", Min, Max);
            plotOne("L2L3res", "BCD", Min, Max);
            outFile = OutFile + "";
            plotOne("L2L3res", "EF", Min, Max);
            plotOne("L2L3res", "GH", Min, Max);

            Min = 0.97, Max = 1.37;
            plotOne("L2rel", "BCD", Min, Max);
            plotOne("L2rel", "EF", Min, Max);
            plotOne("L2rel", "GH", Min, Max);
            plotOne("L2rel", "MC", Min, Max);

            Min = 0.96, Max = 1.04;
            plotOne("L3abs", "BCD", Min, Max);
            plotOne("L3abs", "EF", Min, Max);
            plotOne("L3abs", "GH", Min, Max);
            outFile = OutFile + ")";
            plotOne("L3abs", "MC", Min, Max);


        }

        void PlottSingle(TH2D *hCHSp, TH2D *hCHSn, TH2D *hPUPPIp, TH2D *hPUPPIn, double Min=0.89, double Max=1.2) {
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


                double l = hCHSp->GetXaxis()->GetBinLowEdge(iS);
                double u = hCHSp->GetXaxis()->GetBinUpEdge(iE);
                double c = hCHSp->GetXaxis()->GetBinCenter(iE);

                TH1D *hCHSpProj   =  dynamic_cast<TH1D*>( hCHSp->ProjectionY(SF("%d",rand()), iS, iE) );
                TH1D *hPUPPIpProj =  dynamic_cast<TH1D*>( hPUPPIp->ProjectionY(SF("%d",rand()), iS, iE) );
                TH1D *hCHSnProj   =  dynamic_cast<TH1D*>( hCHSn->ProjectionY(SF("%d",rand()), iS, iE) );
                TH1D *hPUPPInProj =  dynamic_cast<TH1D*>( hPUPPIn->ProjectionY(SF("%d",rand()), iS, iE) );

                hCHSpProj->Draw("hist e ][");
                hCHSpProj->SetLineColor(1);

                hPUPPIpProj->Draw("hist e same ][");
                hPUPPIpProj->SetLineColor(2);

                hCHSnProj->Draw("hist e same ][");
                hCHSnProj->SetLineColor(1);
                hCHSnProj->SetLineWidth(2);
                hCHSnProj->SetLineStyle(2);

                hPUPPInProj->Draw("hist e same ][");
                hPUPPInProj->SetLineColor(2);
                hPUPPInProj->SetLineWidth(2);
                hPUPPInProj->SetLineStyle(2);


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
                    leg->AddEntry(hCHSpProj, "CHS pos", "l");
                    leg->AddEntry(hCHSnProj, "CHS neg", "l");
                    leg->AddEntry(hPUPPIpProj, "PUPPI pos", "l");
                    leg->AddEntry(hPUPPInProj, "PUPPI neg", "l");

                    leg->Draw();
                }

        }
        DrawLatexUp(can->GetPad(1), can->GetPad(7), 1.5, hCHSp->GetTitle() );
    
        can->Print(outFile);
        can->Clear();
        delete can;
    }

};

void JECplotterChsPuppi()
{

    JECPlotter jecPlotter;

    jecPlotter.InitHistos();
    jecPlotter.FillCorr();
    jecPlotter.PlotAll();




}
