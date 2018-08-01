#!/usr/bin/env python
#source /cvmfs/cms.cern.ch/crab3/crab.sh

#File tag for the outputFileNames and jobNames
fileTag = 'puppi0_'

from dataSets  import *

sets =  data16

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    for dataset in sets:
        tagSp = dataset.split('/')
        myName = (fileTag+tagSp[1]+'_'+tagSp[2])[:90]
        #print myName
        #continue

        from WMCore.Configuration import Configuration
        config = Configuration()

        config.section_("User")
        config.User.voGroup = 'dcms'

        config.section_("General")
        config.General.requestName = myName
        config.General.workArea = 'CRAB'

        config.section_("JobType")
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'treeProducerCrab.py'
        config.JobType.pyCfgParams = ['inputFiles='+dataset]
        config.JobType.allowUndistributedCMSSW = True
        config.JobType.maxMemoryMB  = 2500
        config.JobType.maxJobRuntimeMin = 1315

        config.section_("Data")
        config.Data.splitting = 'LumiBased' #LumiBased'
        config.Data.unitsPerJob = 90
        #config.Data.totalUnits  = 10

        #JSON files from https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/
        if '/JetHT/Run2016' in dataset:
            config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' #ReReco
            #config.Data.lumiMask = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'
        elif '/JetHT/Run2017' in dataset:
            config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'

        #config.Data.runRange = '271036-284044'

        config.Data.inputDataset = dataset
        config.Data.outputDatasetTag = myName

        config.section_("Site")
        config.Site.storageSite = "T2_DE_DESY"


        print config.General.requestName, dataset
        crabCommand('submit', config = config)
