#ifndef _utils_
#define _utils_

#include <experimental/filesystem>
#include <vector>
#include <iostream>
#include  <cassert>
#include  "TString.h"

using namespace std::experimental::filesystem;

inline std::vector<TString> GetFiles (const TString source)
{
    std::vector<TString> files;
    path dir(source.Data()); 
    for (auto & file: recursive_directory_iterator(dir)) {
        if (is_regular_file(file.status()) && file.path().extension() == ".root") 
            files.push_back(file.path().c_str());
    }
    return files;
}

inline std::pair<int, int> splitRange(int N, int nSplit, int nNow)
{
    int nStart = N*(nNow+0.)/nSplit; //Included
    int nEnd   = N*(nNow+1.)/nSplit; //notIncluded
    return std::make_pair(nStart, nEnd);
}

inline void PrintCounterAndIncrement(const TString &tag,  const std::pair<int,int> &range)
{
    static long long i =0;
    static long long percentOld = -1;
    static long long nEvSlice = range.second - range.first;

    auto percent = (100ll*i)/nEvSlice;

    if(percent != percentOld && percent % 10 == 0) {
        std::cout << tag <<'\t' << percent << '%' << std::endl;
        percentOld = percent;
    }
    ++i;
}


inline char getPer(int run)
{
     if (273158 <= run && run <= 275376) return 'B';
else if (275657 <= run && run <= 276283) return 'C';
else if (276315 <= run && run <= 276811) return 'D';
else if (276831 <= run && run <= 277420) return 'E';
else if (277981 <= run && run <= 278808) return 'F';
else if (278820 <= run && run <= 280385) return 'G';
else if (281613 <= run && run <= 284044) return 'H';
     else {
         std::cout << "Wrong run " << std::endl;
         std::exit(1);
         return 'x';
     }
}




#endif

