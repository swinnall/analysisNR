import pandas as pd
import csv
import os
import sys
import config


def getFile(path,nSkip,delim,header,names):

    """
    purpose:
        return a pandas dataframe of a given input file

    input:
         path = string of the intended file
         nSkip = int of number of header rows to be skipped
         delim = string of the delimiting separator between values in the input file
         header = indicator of where there is a header present (0) or not (None)
         names = if there is no header, this list of strings is used

    output:
        pandas dataframe

    """
    

    df = pd.read_csv(path, 
                skiprows=nSkip, 
                #sep=delim, # assuming it's automatically delimted now 
                header=header, 
                names=names,
                comment='#',
                na_values =' ',
                skip_blank_lines=True,
                encoding = "utf-8",
                )
    
    # need to replace the dataframe with the whitespace delimeted version
    # as ILL .mft files are not regularly delimted 
    if names != None and 'Q_exp' in names and df.isnull().values.any() == True:

        df = pd.read_csv(path, 
                    skiprows=nSkip, 
                    delim_whitespace=True,
                    header=header, 
                    names=names,
                    comment='#',
                    na_values =' ',
                    skip_blank_lines=True,
                    encoding = "utf-8", 
                    )

    return df
    


def getSampleInfo():

    """
    purpose:
        to get the sample information

    input:
        input paths are defined in config.py by the user under 'paths' heading

    output:
         title = string of the title of the analysis, used in all save outputs
         sampleInfo = pandas dataframe containing the below, one row per sample
             filePaths = dict{key=contrast} of path strings for the raw experimental data
             inputLabel = dict{key=contrast} of label strings for the raw experimental data
             contrastList = list of strings defining each contrast
             lipids = dict{key=contrast:string} of component lipids
             ratios = dict{key=contrast:string} of ratios of the component lipids
    """

    # get the title of sample info file
    with open(config.sample_info_fname, newline='') as f:
        title = list(csv.reader(f))[0][0].split('=')[1]


    # get the sample info file
    sampleInfo = getFile(path=config.sample_info_fname, nSkip=1, delim=',',header=0,names=None)


    # check that a file has been inputted
    if len(sampleInfo) == 0:
        print(f'\nFatal Error: No input files in {config.sample_info_fname}')
        print('Possible Reason: All rows are commented out.')
        sys.exit()


    # for each file get the full path
    file_path_list = []
    for fname in sampleInfo['fname']:

        # initialise string, only gets updated if file is found
        file_path_str = ''

        # find full path of sample data and store in filePaths
        for root,dirs,files in os.walk('../input/'):
            
            # iterate through all of the file names
            for name in files:

                # if a match exists then store
                if name in [fname+".txt",fname+".dat",fname+".mft"]:
                    file_path_str = os.path.join(root,name)
                    file_path_list.append(file_path_str)


        # if file path has not been updated then the file has not been found
        if file_path_str == '':
            print('\nFatal Error: file not found.')
            print('Possible Reason: typo or non txt/dat file extension.')
            print(f'File name = {fname}')
            sys.exit()

    if len(file_path_list) > len(sampleInfo):
        print(f'\nFatal Error: File(s) found ({len(file_path_list)}) exceed described samples ({len(sampleInfo)}).')
        print('List of files found:')
        print(*file_path_list,sep='\n')
        sys.exit()

    if len(file_path_list) < len(sampleInfo):
        print(f'\nFatal Error: File(s) found ({len(file_path_list)}) is less than defined samples ({len(sampleInfo)}).')
        print('List of files found:')
        print(*file_path_list,sep='\n')
        sys.exit()

    # load into new df column
    sampleInfo['filePath'] = file_path_list

    # update to uppercase
    sampleInfo['contrast'] = sampleInfo['contrast'].str.upper()
    sampleInfo['lipid'] = sampleInfo['lipid'].str.upper()

    return title, sampleInfo



def getExperimentalData(sampleInfo):

    """
    purpose:
        to read the sample information dataframe and store the raw experimental
        data in the sampleInfo dataframe
    """

    # initialise variables
    Q_exp_list, R_exp_list, R_err_exp_list, qmin_list, qmax_list = [], [], [], [], []


    # for each file get the experimental data
    for path_ in sampleInfo['filePath']:

        # get the file extension (last element in split string)
        suffix = path_.split('.')[-1]

        # select the correct delimeter
        if suffix == 'txt':
            delim_ = '\t'
            names_ = ['Q_exp','R_exp','R_err_exp','Q_exp_err']

        elif suffix == 'dat':
            delim_ = ' '
            names_ = ['Q_exp','R_exp','R_err_exp']

        elif suffix == 'mft':
            delim_ = ' '
            names_ = ['Q_exp','R_exp','R_err_exp', 'Q_exp_err']


        # get the experimental data as a data frame
        experiment_df = getFile(path=path_, nSkip=1, delim=delim_, header=None,\
                                names=names_)

        #print(experiment_df); sys.exit()
        
        # truncate experiment_df to only include low Q data
        if config.surface_excess_mode == True:
            experiment_df = experiment_df[experiment_df['Q_exp'] < 0.07]


        # store each series in the list
        Q_exp_list.append(experiment_df['Q_exp'])
        R_exp_list.append(experiment_df['R_exp'])
        R_err_exp_list.append(experiment_df['R_err_exp'])
        qmin_list.append(min(experiment_df['Q_exp']))
        qmax_list.append(max(experiment_df['Q_exp']))


    # load into new df column
    sampleInfo['Q_exp'] = Q_exp_list
    sampleInfo['R_exp'] = R_exp_list
    sampleInfo['R_err_exp'] = R_err_exp_list

    sampleInfo['qmin'] = qmin_list
    sampleInfo['qmax'] = qmax_list

    return sampleInfo
