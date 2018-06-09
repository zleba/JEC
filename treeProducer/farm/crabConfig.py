#source /cvmfs/cms.cern.ch/crab3/crab.sh

from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'JECa_RunB'
config.General.workArea = 'radekTestE'

config.section_("JobType")
config.JobType.pluginName = 'ANALYSIS'
config.JobType.psetName = 'treeProducerNew.py'
#config.JobType.psetName = 'flatData-TTJets-FullTriggers_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB  = 2500
config.JobType.maxJobRuntimeMin = 1315

config.section_("Data")
#config.Data.inputDataset = '/JetHT/Run2016G-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-03Feb2017-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2016B-03Feb2017_ver2-v2/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2016B-07Aug17_ver2-v1/MINIAOD'

#config.Data.inputDataset = '/JetHT/Run2016D-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016E-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-03Feb2017-v1/MINIAOD'
config.Data.splitting = 'LumiBased' #LumiBased' 
config.Data.unitsPerJob = 20
#config.Data.totalUnits  = 10
config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
#config.Data.runRange = '271036-284044' 
config.Data.outputDatasetTag = 'JECa_RUNB'


config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"

#if __name__ == '__main__':
#    from CRABAPI.RawCommand import crabCommand
#    crabCommand('submit', config = config)


sets = [
'/JetHT/Run2016B-07Aug17_ver2-v1/MINIAOD',
'/JetHT/Run2016C-07Aug17-v1/MINIAOD',
'/JetHT/Run2016D-07Aug17-v1/MINIAOD',
'/JetHT/Run2016E-07Aug17-v1/MINIAOD',
'/JetHT/Run2016F-07Aug17-v1/MINIAOD',
'/JetHT/Run2016G-07Aug17-v1/MINIAOD',
'/JetHT/Run2016H-07Aug17-v1/MINIAOD',
]

sets = [ sets[MYID] ]

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    for dataset in sets:
        config.Data.inputDataset = dataset
        n = 'JECc_'+dataset.split('/')[2]
        config.Data.outputDatasetTag = n
        config.General.requestName   = n
        print config.General.requestName
        crabCommand('submit', config = config)
