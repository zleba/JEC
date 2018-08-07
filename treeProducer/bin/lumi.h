#ifndef lumi_h
#define lumi_h

#include <map>
#include <vector>
#include <fstream>
#include "TString.h"

static TString lumiPath = "/nfs/dust/cms/user/connorpa/SMPJ/effective_luminosity_Run2016BtoH/";

class Luminosity {
    vector<map<int, double>> lumiMap;


    public:

    static map<int,double> GetLumis(TString fName)
    {
        std::ifstream file(fName);
        assert(file.good());
        
        map<int,double> lumiMap;
        while(1) {
            TString trigName;
            double lum;
            file >> trigName >> lum;
            if(!file.good()) break;
            trigName.ReplaceAll("HLT_PFJet", "");
            trigName = trigName(0, trigName.Length()-3);

            int trigTr = trigName.Atoi();
            if(trigTr != 500) {
                lumiMap[trigTr] += lum;
            }
        }
        file.close();
        return lumiMap;

    }

    void LoadLumis() {
        vector<vector<TString>> files;
        files.push_back({"processedLumis2016B.json.txt"});
        files.push_back({"processedLumis2016C.json.txt"});
        files.push_back({"processedLumis2016D.json.txt"});
        files.push_back({"processedLumis2016E.json.txt"});
        files.push_back({"processedLumis2016Fearly.json.txt", "processedLumis2016Flate.json.txt" }); //Merging F into one period
        files.push_back({"processedLumis2016G.json.txt"});
        files.push_back({"processedLumis2016H.json.txt"});

        lumiMap.resize(files.size() + 1); //Index 0 is dedicated for x-section summed over all periods


        //Read luminosities to map
        for(unsigned i = 0; i < files.size(); ++i) {
            for(unsigned j = 0; j < files[i].size(); ++j) {
                auto lum = GetLumis(lumiPath + files[i][j]);
                for(auto v : lum) {
                    lumiMap[i+1][v.first] += v.second;
                    //cout << "Printout " << v.first <<" "<< v.second << endl;
                }
            }
        }

        //Evaluate the total luminosity
        for(unsigned i = 1; i < lumiMap.size(); ++i)
            for(auto v : lumiMap[i])
                lumiMap[0][v.first] += v.second;

        /*
        //Sum B,C,D and make them identical
        map<int,double> lumiMapSingle;
        for(auto v : lumiMap[0]) {
            double sum = lumiMap[1][v.first] + lumiMap[2][v.first] + lumiMap[3][v.first];
            lumiMap[1][v.first] = lumiMap[2][v.first] = lumiMap[3][v.first] = sum;
        }
        */



        //Assert identical trigger configuration:
        for(unsigned i = 1; i < lumiMap.size(); ++i) {
            for(auto m : lumiMap[0])
                assert(lumiMap[i].count(m.first));
            //assert(lumiMap[0].size() == lumiMap[i].size());
        }


        for(auto v : lumiMap[0]) {
            cout << v.first <<" "<< v.second << endl;
        }


    }

    //runLum, totLum, trigID
    tuple<double,double,int> GetWeightID(int periodId, double pT0)
    {
        //Trigger tresholds from Patrick
        static const map<int,int> trigTrsh = {
            {40 , 74},
            {60 , 84},
            {80 , 114},
            {140 , 196},
            {200 , 245},
            {260 , 330},
            {320 , 395},
            {400 , 468},
            {450 , 507},
        };

        assert(periodId < int(lumiMap.size()));
        int id = trigTrsh.size() - 1;
        for(auto it = trigTrsh.rbegin(); it != trigTrsh.rend(); ++it, --id) 
            if(pT0 >= it->second) {
                return make_tuple(1. / lumiMap[0].at(it->first),  1. / lumiMap[periodId].at(it->first), id);
            }

        return make_tuple(0.,0., -1);

    }

};

#endif
