
data16= [
        '/JetHT/Run2016B-07Aug17_ver2-v1/MINIAOD',
        '/JetHT/Run2016C-07Aug17-v1/MINIAOD',
        '/JetHT/Run2016D-07Aug17-v1/MINIAOD',
        '/JetHT/Run2016E-07Aug17-v1/MINIAOD',
        '/JetHT/Run2016F-07Aug17-v1/MINIAOD',
        '/JetHT/Run2016G-07Aug17-v1/MINIAOD',
        '/JetHT/Run2016H-07Aug17-v1/MINIAOD',
        ]

data17= [
        '/JetHT/Run2017B-31Mar2018-v1/MINIAOD',
        '/JetHT/Run2017C-31Mar2018-v1/MINIAOD',
        '/JetHT/Run2017D-31Mar2018-v1/MINIAOD',
        '/JetHT/Run2017E-31Mar2018-v1/MINIAOD',
        '/JetHT/Run2017F-31Mar2018-v1/MINIAOD',
        ]

import os,re

#Get the short name
def Short(s):
    return re.search(r'[^/]*(Run|QCD)[^/]*/[^-]*', s).group()


#addresses of all files in give set (by default all are given)
def GetFiles(s, limit=0):
    return os.popen('dasgoclient --query="file dataset='+s+'" -limit='+str(limit)).read().splitlines()


