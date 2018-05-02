#ifndef _QCDjet_
#define _QCDjet_

#include <vector>
//#include "TObject.h"
#include "Math/PtEtaPhiM4D.h"

using namespace std;

class QCDjet {
    public:
    int flavor      ;
    int flavorHadron;
    float chf       ;
    float nhf       ;
    float phf       ;
    float elf       ;
    float muf       ;
    int chm         ;
    int nhm         ;
    int phm         ;
    int elm         ;
    int mum         ;
    float jetJECtot;
    float jetJECl2l3Res;

    float  btag;
    float area;
    float unc;
    ROOT::Math::PtEtaPhiM4D<float> p4;
    bool tightID;

    QCDjet() {
        /*
        flavor      = 0;
        flavorHadron= 0;
        chf       = 0;
        nhf       = 0;
        phf       = 0;
        elf       = 0;
        muf       = 0;
        chm         = 0;
        nhm         = 0;
        phm         = 0;
        elm         = 0;
        mum         = 0;
        jetJECfact= 0;

        isBtag= 0;
        btag= 0;
        area= 0;
        unc= 0;
        p4 = TLorentzVector();
        tightID= 0;
        */
    }
    virtual ~QCDjet(){}

    //ClassDef(QCDjet,1)

};


#if defined(__ROOTCLING__)
#pragma link C++ class QCDjet +;
#pragma link C++ class vector<QCDjet> +;

#pragma link C++ class HLTjet +;
#pragma link C++ class vector<HLTjet> +;
#pragma link C++ class vector<ROOT::Math::PtEtaPhiM4D<float> > +;
#endif

#endif
