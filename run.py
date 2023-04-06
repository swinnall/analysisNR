"""
purpose: 
    run script for fitting experimental neutron data 
    
input: 
    sample-info.txt for experimental file name corresponding sample info
    config.py for parameter information
    raw data files from ../input 

analysis:
    subscript_fitNR module processes the inputs and generates model data 

output: 
    subscript_genFig plots experimental and model data 
    writeParams function stores global objective parameters in a txt file
    varied and constrained parameters are also printed to the terminal 
    
to do: 
    update figure to full class representation 
"""

# ---------------------------------------------------------------------------#
# IMPORTS

import config
from subscript_getFiles import getSampleInfo, getExperimentalData
from subscript_fitNR import fitModel
from subscript_genFig import Figure


# ---------------------------------------------------------------------------#
# READ SAMPLE INFORMATION

# get sample instructions
title, file_paths, inputLabels, contrastList, lipids, ratios = getSampleInfo()


# ---------------------------------------------------------------------------#
# GET EXPERIMENTAL DATA 
Q, expNR, expNR_err, labels, qmin, qmax = getExperimentalData(file_paths, \
                        inputLabels, contrastList, nContrasts=len(contrastList))

    
# ---------------------------------------------------------------------------#
# GENERATE MODEL DATA 
    
# generate model data
modelQ, modelNR, global_objective, reduced_chisq, APM = fitModel(title, \
                          file_paths, contrastList, qmin, qmax, lipids, ratios)

    
# ---------------------------------------------------------------------------#
# OUTPUT 

fig = Figure()

if config.surface_excess_mode == True: 
            
    fig.plotExcess(contrastList, global_objective)
    
    
else:
            
    fig.plotStruct(Q,expNR,expNR_err,modelQ,modelNR,labels,title,row=1,col=1,showlegend=False)    
    

    print('\n\n--------------------------------------------------------------------')

    print('\nVaried Parameters:')
    for par_vary in global_objective.parameters.varying_parameters():
        print(par_vary.name, "{:.1f}".format(par_vary.value), par_vary.bounds)
    
    print('\nConstrained Parameters:')
    for par_constrain in global_objective.parameters.constrained_parameters():
        print(par_constrain.name, "{:.1f}".format(par_constrain.value))
      
        
    # print calculated quantities
    print('\nReduced chisqr = %.1f' %reduced_chisq)
    print(f'\nArea per molecule = {APM:.1f} A^2')


fig.save_figure(title)

fig.show_figure()
