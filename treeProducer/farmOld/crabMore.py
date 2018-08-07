sets = [
'/JetHT/Run2016B/JetHT/MINIAOD/07Aug17_ver2-v1',
'/JetHT/Run2016C/JetHT/MINIAOD/07Aug17-v1',
'/JetHT/Run2016D/JetHT/MINIAOD/07Aug17-v1',
'/JetHT/Run2016E/JetHT/MINIAOD/07Aug17-v1',
'/JetHT/Run2016F/JetHT/MINIAOD/07Aug17-v1',
'/JetHT/Run2016G/JetHT/MINIAOD/07Aug17-v1',
'/JetHT/Run2016H/JetHT/MINIAOD/07Aug17-v1',
]


for dataset in sets:
    #config.Data.inputDataset = dataset
    #config.Data.outputDatasetTag = config.General.requestName = 'JEC_'+dataset.split('/')[2]
    name = 'JEC_'+dataset.split('/')[2]
    #print config.General.requestName
    from subprocess import call
    call(["echo", 'submit', 'Data.inputDataset='+ dataset, 'General.requestName='+name,  'Data.outputDatasetTag='+name])

