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
title, sampleInfo = getSampleInfo()


# ---------------------------------------------------------------------------#
# GET EXPERIMENTAL DATA
sampleInfo = getExperimentalData(sampleInfo)


# ---------------------------------------------------------------------------#
# GENERATE MODEL DATA

# generate model data
sampleInfo, global_objective, reduced_chisq, APM, sigma = fitModel(title,sampleInfo)


# ---------------------------------------------------------------------------#
# OUTPUT

fig = Figure()

if config.surface_excess_mode == True:

    fig.plotExcess(sampleInfo, global_objective)


else:

    fig.plotStruct(title,sampleInfo,row=1,col=1,showlegend=False)


    print('\n\n--------------------------------------------------------------------')

    # print reminders for when merging third and headgroup layers
    if config.with_drug_layer2 == True and config.d3_vary == False:
        print('\nReminder: user has fixed d3 while merging with headgroup layer.')
        print('This introduces the constraint: d3 = config.d3_0 - d2')

    if config.with_drug_layer2 == True and config.d3_vary == True:
        print('\nReminder: user has fitted d3 while merging with headgroup layer.')
        print('This means the total drug thickness = d3 + d2')



    print('\nVaried Parameters:')
    for par_vary in global_objective.parameters.varying_parameters():

        if par_vary.name == 'bkg':
            print(par_vary.name, "{:.2e}".format(par_vary.value), par_vary.bounds)
        else:
            print(par_vary.name, "{:.2f}".format(par_vary.value), par_vary.bounds)


    print('\nConstrained Parameters:')
    for par_constrain in global_objective.parameters.constrained_parameters():
        print(par_constrain.name, "{:.2f}".format(par_constrain.value))


    # print calculated quantities
    print('\nReduced chisqr = %.2f' %reduced_chisq)
    print(f'\nArea per molecule = {APM:.2f} A^2')
    print(f'\nCharge density = {sigma:.2f} C/m^2\n NB: Multiply by % charged for correct value')



fig.save_figure(title)

fig.show_figure()
