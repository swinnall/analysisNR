import pandas as pd
import csv
import config


def getFile(path,nSkip,delim):
    
    """
    purpose: 
        return a pandas dataframe of a given input file 
        
    input: 
         path = string of the intended file 
         nSkip = int of number of header rows to be skipped
         delim = string of the delimiting separator between values in the input file

    output: 
        pandas dataframe 
         
    """  
    
    return pd.read_csv(path, skiprows=nSkip, sep=delim, comment='#', \
                       na_values =' ', skip_blank_lines=True, \
                       encoding = "utf-8") 


def getSampleInfo():

    """
    purpose: 
        to get the sample information 
        
    input: 
        input paths are defined in config.py by the user under 'paths' heading

    output: 
         title = string of the title of the analysis, used in all save outputs 
         filePaths = dict{key=contrast} of path strings for the raw experimental data 
         inputLabel = dict{key=contrast} of label strings for the raw experimental data
         contrastList = list of strings defining each contrast
         lipids = dict{key=contrast:string} of component lipids
         ratios = dict{key=contrast:string} of ratios of the component lipids
    """     

    # variable definitions 
    filePaths    = {}
    inputLabels  = {}
    lipids       = {}
    ratios       = {}
    contrastList = []
    
    # get the title of sample info file 
    with open(config.sample_info_str, newline='') as f:
        title = list(csv.reader(f))[0][0].split('=')[1]
    
    # get the sample info file 
    sampleInfoFile = getFile(path=config.sample_info_str, nSkip=1, delim=',')
    
    # calculate number of files in sample info
    nFiles = len(sampleInfoFile)
    
    # for each file get the sample specific information and store in dicts
    for i in range(nFiles):

        contrast = sampleInfoFile["contrast"][i]
        contrastList.append(contrast)

        filePaths[contrast] = "../input/" + config.input_sub_str + '/' + sampleInfoFile["fname"][i] + ".txt"
        inputLabels[contrast] = sampleInfoFile["label"][i]

        lipids[contrast] = sampleInfoFile["lipid"][i]
        ratios[contrast] = sampleInfoFile["ratio"][i]

    return title, filePaths, inputLabels, contrastList, lipids, ratios



def getExperimentalData(file_paths, inputLabels, contrastList, nContrasts):
    
    """
    purpose: 
        to read the sample information dataframe and store the raw experimental
        data 
        
    input: 
         filePaths = dict{key=contrast} of path strings for the raw experimental data 
         inputLabel = dict{key=contrast} of label strings for the raw experimental data
         contrastList = list of strings defining each contrast
         nContrasts = int of the total number of contrasts 

    output: 
         dictionaries of the experimental quantities accessed via contrast key
    """  
    
    # initialise variables
    expQ = {}
    expNR = {}
    expNR_err = {}
    labels = {}
    qmin = {}
    qmax = {}
    
    # for each contrast, store the experimental information
    for i in range(nContrasts):
        contrast = contrastList[i]

        if file_paths.get(contrast) is not None:
            experimentFile = getFile(path=file_paths.get(contrast), nSkip=1, delim='\t')

            expQ[i] = experimentFile[experimentFile.columns.values[0]]
            expNR[i] = experimentFile[experimentFile.columns.values[1]]
            expNR_err[i] = experimentFile[experimentFile.columns.values[2]]
            labels[i] = inputLabels.get(contrast)
            
            qmin[contrast] = min(expQ.get(i))
            qmax[contrast] = max(expQ.get(i))
    
    return expQ, expNR, expNR_err, labels, qmin, qmax

