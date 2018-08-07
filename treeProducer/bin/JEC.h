#ifndef JECs_h
#define JECs_h

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "JEC/treeProducer/interface/QCDjet.h"

using namespace edm;
using namespace std;

static const string JECpath = "/nfs/dust/cms/user/connorpa/DAS18/JEStables";


class JECs  {

    bool isMC;
    FactorizedJetCorrector *jecL1Fast, *jecL2Relative, *jecL3Absolute, *jecL2L3Residual;

    double ApplyJECL1Fast (double &pt, double eta, double area, double rho)
    {
        // ---- Evaluating the L1Fast correction factor ---- //
        jecL1Fast->setJetPt(pt);
        jecL1Fast->setJetEta(eta);
        jecL1Fast->setJetA(area);
        jecL1Fast->setRho(rho);
        double corr = jecL1Fast->getCorrection();
        pt *= corr;
        return corr;
    }

    double ApplyJECL2Relative (double &pt, double eta)
    {
        // ---- Evaluating the L2Relative correction factor ---- //
        jecL2Relative->setJetPt(pt);
        jecL2Relative->setJetEta(eta);
        double corr = jecL2Relative->getCorrection();
        pt *= corr;
        return corr;
    }

    double ApplyJECL3Absolute (double &pt, double eta)
    {
        // ---- Evaluating the L2Relative correction factor ---- //
        jecL3Absolute->setJetPt(pt);
        jecL3Absolute->setJetEta(eta);
        double corr = jecL3Absolute->getCorrection();
        pt *= corr;
        return  corr;
    }

    double ApplyJECL2L3Residual (double &pt, double eta)
    {
        jecL2L3Residual->setJetPt(pt);
        jecL2L3Residual->setJetEta(eta);
        double corr =  jecL2L3Residual->getCorrection();
        pt *= corr;
        return corr;
    }


public:

    double getJEC(const QCDjet & jet, double  rho)
    {
        double pt = jet.p4.Pt();
        double eta = jet.p4.Eta();
        double area = jet.area;

        double corr = 1;
        //Apply JEC
        corr *= ApplyJECL1Fast(pt, eta, area, rho);
        corr *= ApplyJECL2Relative(pt, eta);
        corr *= ApplyJECL3Absolute(pt, eta);

        if(!isMC)
            corr *= ApplyJECL2L3Residual(pt, eta);

        return corr;
    }

    JECs (string name, string jetType) :
            isMC(name.find("MC") != string::npos),
            jecL1Fast      (nullptr),
            jecL2Relative  (nullptr),
            jecL3Absolute  (nullptr),
            jecL2L3Residual(nullptr)
    {
//AK4PFchs
        auto initCorr = [&](string corrType) {
            string wholePath =JECpath + "/" + name + "/" + name + "_" + corrType +"_"+jetType+".txt";
            //cout << "WholePath " << wholePath  << endl;
            auto corPars       = new JetCorrectorParameters(wholePath);
            vector<JetCorrectorParameters> corParsVec;
            corParsVec    .push_back(*corPars);
            return (new FactorizedJetCorrector(corParsVec));
        };

        jecL1Fast = initCorr("L1FastJet");
        jecL2Relative = initCorr("L2Relative");
        jecL3Absolute = initCorr("L3Absolute");
        if(!isMC) jecL2L3Residual = initCorr("L2L3Residual");

    }

};


#endif
