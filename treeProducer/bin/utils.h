#ifndef _utils_
#define _utils_

#include <experimental/filesystem>
#include <vector>

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

#endif

