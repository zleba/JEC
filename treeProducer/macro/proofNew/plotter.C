#define SF TString::Format

#include "/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/libPlottingHelper.so)

#include "TH3D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TStyle.h"



using namespace PlottingHelper;

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
    vector<TH2D *> hEtaPtPUPPI;
    vector<TH2D *> hEtaPtPUPPIalone;



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
        hEtaPtPUPPI.resize(nPer);
        hEtaPtPUPPIalone.resize(nPer);


        for(int i = 0; i < nPer; ++i) {
            char per = 'A' + i;
            for(int j = 0; j < sizeof(hBalEtaPtAll)/sizeof(hBalEtaPtAll[0]); ++j) {
                hBalEtaPtAll[j][i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt%d_%c", j, per)));
                //hBalEtaPtAll[j][i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt_%c",  per)));
                assert(hBalEtaPtAll[j][i]);
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

    void MeanAsym(int shift=0, TString style="");
    void JEC();
    void PtEtaDep();
    void Unmatched();
    void MatchingFactorsTimeDep(int rho=0);
    void MatchingFactorsPileUpDep();

};

pair<TPad*, TPad*> TitleSpace(int perID)
{
    TVirtualPad *can = gPad;
    const double sp = 0.97;
    TPad *up = new TPad(SF("%d",rand()),"",   0, sp, 1, 1.0);

    up->Draw();

    up->cd();
    TLatex *lat = new TLatex;
    TString s;
    if(perID != 0) s =  SF("Run 2016%c", perID + 'A');
    else s =  SF("Run 2016 all");
    lat->SetTextSize(0.6);
    lat->DrawLatexNDC(0.8, 0.3, s);


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

    cout << "I am here " << __LINE__ << endl;
    dnPad->Divide(5, 5, 0, 0);
    dnPad->SetBottomMargin(0.1);
    cout << "I am here " << __LINE__ << endl;
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

    cout << "I am here " << __LINE__ << endl;


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

void PLOTTER::Unmatched()
{
    auto hPUPPI = hEtaPtPUPPI[perID];
    auto hPUPPIa = hEtaPtPUPPIalone[perID];

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

        //hChs->Scale(1., "width");

        //hPuppiA->SetMinimum(0.00);
        //hPuppiA->SetMaximum(0.20);

        hPuppiA->Draw();
        //hPuppi->Draw("same");
        //hChs->Draw("same");

        //SetFontSizes(hPuppiA, 0.07);

        GetXaxis()->SetRangeUser(34, 1780);
        GetYaxis()->SetRangeUser(0.0, 0.029);
        GetXaxis()->SetMoreLogLabels();
        GetXaxis()->SetNoExponent();

        SetFTO({12}, {5.1}, {1.63, 3.0, 0.3, 5.4});

        GetFrame()->SetTitle("");
        if(i == 1)
            GetYaxis()->SetTitle("#sigma^{unmatched}_{PUPI}/#sigma_{PUPPI}");
        if(i == 20)
            GetXaxis()->SetTitle("p_{T}^{PUPPI} [GeV]");

        DrawLatexUp(-1, SF("%1.1f < |#eta| < %1.1f", l, u));
    }
    
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

void PLOTTER::MatchingFactorsPileUpDep()
{

    TH3D *hBalEtaPt;
    hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAll[0][0]->Clone(SF("%d",rand())));
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

        vector<TProfile*> prof(5);
        for(int k = 0; k < 5; ++k) {
            iS = min(19,iS);
            iE = min(19,iE);
            hBalEtaPtAll[k][0]->GetXaxis()->SetRange(iS,iE);
            TH2D *hTemp = dynamic_cast<TH2D*>(hBalEtaPtAll[k][0]->Project3D(SF("%d_yz",rand()))); 
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
            //for(int k = 0; k < 5; ++k) {
                //TString t = k != 0 ? SF("Run %c",'A'+ k) : "All Runs";
                //leg->AddEntry(prof[k], t, "l");
            //}
            leg->AddEntry(prof[0], "any #rho", "l");
            leg->AddEntry(prof[1], "#rho < 15", "l");
            leg->AddEntry(prof[2], "15 < #rho < 22", "l");
            leg->AddEntry(prof[3], "22 < #rho < 30", "l");
            leg->AddEntry(prof[4], "30 < #rho", "l");
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












void PLOTTER::MeanAsym(int shift, TString style)
{

    TH3D *hBalEtaPt;
    if(shift == 0)
        hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAll[0][perID]->Clone(SF("%d",rand())));
    else if(shift == 1)
        hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAllUp[perID]->Clone(SF("%d",rand())));
    else
        hBalEtaPt = dynamic_cast<TH3D*>( hBalEtaPtAllDn[perID]->Clone(SF("%d",rand())));

    assert(hBalEtaPt);

    hBalEtaPt->RebinY(4);
    gStyle->SetPaintTextFormat("2.2f");
    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);


    dnPad->SetMargin(0.1, 0.1, 0.1, 0.1);
    dnPad->SetLogx();

    TProfile2D * hProf = hBalEtaPt->Project3DProfile();
    hProf->GetXaxis()->SetRangeUser(20, 2000);
    hProf->BuildOptions(0.8, 1.2, "s");

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


    if(style == "colz")
        hProf->Draw("colz");
    else
        hProf->Draw("text");

    can->Print(outName);
    can->Clear();
    delete hProf;
    delete can;
}







void plotter()
{
    PLOTTER plot;
    //TString myOut = "plots/Spring16_25nsV6.pdf";
    //plot.Init("histos/Spring16_25nsV6__Spring16_25nsV6.root", myOut+"(");

    TString myOut = "plots/Summer16_07Aug2017V5newnew.pdf";
    plot.Init("histos/Summer16_07Aug2017V5__Summer16_07Aug2017V5newnew.root", myOut+"(");
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

    plot.perID = 0;
    plot.Unmatched();
    plot.outName = myOut;
    plot.AsymmetryEtaPtTimeDep(1, 5);
    plot.AsymmetryEtaPtTimeDep(6, 10);
    plot.AsymmetryEtaPtTimeDep(11, 15);
    plot.AsymmetryEtaPtTimeDep(16, 18);
    plot.MatchingFactorsTimeDep(0);
    //plot.MatchingFactorsTimeDep(1);
    //plot.MatchingFactorsTimeDep(2);
    //plot.MatchingFactorsTimeDep(3);


    plot.AsymmetryEtaPtPileUpDep(1, 5);
    plot.AsymmetryEtaPtPileUpDep(6, 10);
    plot.AsymmetryEtaPtPileUpDep(11, 15);
    plot.AsymmetryEtaPtPileUpDep(16, 18);



    plot.outName = myOut + ")";
    plot.MatchingFactorsPileUpDep();

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
