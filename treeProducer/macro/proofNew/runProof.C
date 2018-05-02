void AddPythia(TChain *ch);
void AddMadgraph(TChain *ch);
#include <string>
#include <cmath>

#define SF TString::Format 

using namespace std;

const TString ntuplePath = "/nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged";

void runProof(int nMax=1, int nNow=0, char per = 'A') { 

    
    //gSystem->Load("../../plugins/QCDjet.h");

    //gROOT->ProcessLine(".L  /afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/plugins/QCDjet.h");


    TChain * chain = new TChain("ak4/events");
    //AddPythia(chain);
    //AddMadgraph(chain);
    //return;

    if(per == 'A')
        chain->Add(ntuplePath +  "/jets[B-H].root"); 
    else {
        chain->Add(SF(ntuplePath + "/jets%c.root", per)); 
    }


    //bool on = nworkers>0;

    //chain->LoadTree(-1);
    int N = chain->GetEntries();

    if (0) { 
        TProof::Open("workers=11");//"workers=3");
        //TString connect = gSystem->GetFromPipe("pod-info -c");
        //TProof::Open(connect);//"workers=3");
        //
        //TProof::Open(TString::Format("workers=%d",nworkers) );

        chain->SetProof(); 
    }


    cout << "POSITION " << nMax << " "<< nNow << " | " << N << endl;


    int nStart = lround(nNow / double(nMax) * double(N));
    int nEnd   = lround((nNow+1) / double(nMax) * double(N));


    //TStopwatch w; w.Start();
    chain->Process("matching.C+", "", nEnd - nStart + 1, nStart);
    //chain->Process("matching.C+","", 10000, 0);
    //w.Print();

}

