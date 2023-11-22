from refnx.dataset import ReflectDataset
from refnx.analysis import Objective, Parameter
from refnx.reflect import SLD, ReflectModel
from refnx.analysis import CurveFitter, GlobalObjective
from scipy.optimize import NonlinearConstraint
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import integrate
import numpy as np
import functools
import sys
import pandas as pd

from subscript_calcSLD import calc_sld
import config


class Model:

    """
    purpose:
        Model is a class where each method is a step in the refnx procedure
        for generating a structural slab model of the lipid monolayer and drug
        layer where required

    input:
        parameters are defined in config.py by the user

    output:
        an objective (initial parameter solution) is created for each sample
        contrast that is loaded into the class, this is stored for a later
        global fit
    """

    def __init__(self):

        # define roughness
        self.roughness = Parameter(config.rough_0, name='Roughness')

        # define roughness backing
        if config.rough_bkg_incl == True:
            self.roughness_backing = Parameter(config.rough_bkg_0, bounds=(config.rough_bkg_lb,config.rough_bkg_ub), vary=config.rough_bkg_vary, name="Rough_bkg")

        # define front slab (always air)
        self.air = SLD(0.0, name='air')

        # define thicknesses
        self.thickness_tails = Parameter(config.d1_0, bounds=(config.d1_lb,config.d1_ub), vary=config.d1_vary, name='$d_1 (\AA)$')
        self.thickness_heads = Parameter(config.d2_0, bounds=(config.d2_lb,config.d2_ub), vary=config.d2_vary, name='$d_2 (\AA)$')
        self.thickness_layer_3 = Parameter(config.d3_0, bounds=(config.d3_lb,config.d3_ub), vary=config.d3_vary, name='Layer 3 Thick.')
        self.thickness_layer_4 = Parameter(config.d4_0, bounds=(config.d4_lb,config.d4_ub), vary=config.d4_vary, name='Layer 4 Thick.')
        self.thickness_layer_5 = Parameter(config.d5_0, bounds=(config.d5_lb,config.d5_ub), vary=config.d5_vary, name='Layer 5 Thick.')
        self.thickness_layer_6 = Parameter(config.d6_0, bounds=(config.d6_lb,config.d6_ub), vary=config.d6_vary, name='Layer 6 Thick.')
        self.thickness_layer_7 = Parameter(config.d7_0, bounds=(config.d7_lb,config.d7_ub), vary=config.d7_vary, name='Layer 7 Thick.')

        # define linked thickness
        if config.link_drug_thick == True:
            self.linked_drug_thick = Parameter(config.d3_link_0, bounds=(config.d3_link_lb,config.d3_link_ub), vary=config.d3_link_vary, name='Layer 3 Linked Thick.')

        # define tail solvent (always 0, no solvent in tails)
        self.tail_solvent = Parameter(0.0, name='Tail Solv.')

        # define how much solvent in drug layer (constant between contrasts)
        if config.with_drug_layer3 == True:
            self.layer3_solvent = Parameter(config.solv3_0, bounds=(config.solv3_lb,config.solv3_ub), vary=config.solv3_vary, name='Layer 3 Solv.')

        if config.with_drug_layer4 == True:
            self.layer4_solvent = Parameter(config.solv4_0, bounds=(config.solv4_lb,config.solv4_ub), vary=config.solv4_vary, name='Layer 4 Solv.')

        if config.with_drug_layer5 == True:
            self.layer5_solvent = Parameter(config.solv5_0, bounds=(config.solv5_lb,config.solv5_ub), vary=config.solv5_vary, name='Layer 5 Solv.')

        if config.with_drug_layer6 == True:
            self.layer6_solvent = Parameter(config.solv6_0, bounds=(config.solv6_lb,config.solv6_ub), vary=config.solv6_vary, name='Layer 6 Solv.')

        if config.with_drug_layer7 == True:
            self.layer7_solvent = Parameter(config.solv7_0, bounds=(config.solv7_lb,config.solv7_ub), vary=config.solv7_vary, name='Layer 7 Solv.')

        # SLD constant (based on reasonable guess) and solvent fixed at 0
        if config.surface_excess_mode == True:
            self.rho_slab = SLD(config.sld_slab_0, name='SLD_slab_layer')
            self.slab_solvent = Parameter(0.0, name='Solv_slab_layer')


        # define lists for data
        self.exp_list = []
        self.obj_list = []
        self.struct_list = []
        self.model_list = []




    def load_contrast(self, file, contrast, lipid, ratio, title, label, colour):

        # define file / sample information parameters
        self.file = file
        self.contrast = contrast.split('_')

        # receiving a specific lipid and ratio object from main() via lipids.get(i) for e.g.
        self.lipid = lipid
        self.ratio = ratio

        # define the title and label of the figure and experimental data for saving SLD plots
        self.title = title
        self.label = label

        # define the colour for each model
        self.colour = colour

        # extract system information from contrast
        self.subphase = self.contrast[0]




    def make_backing_slab(self):

        if self.subphase == 'ACMW':

            # define initial variable
            subphase_init = config.acmw_0

            # define bounds for varying parameters
            subphase_lb, subphase_ub = config.acmw_lb, config.acmw_ub

            # define boolean for whether parameter is fixed or varying
            subphase_vary = True if config.acmw_vary == True else False


        elif self.subphase == 'D2O':

            # define initial variable
            subphase_init = config.d2o_0

            # define bounds for varying parameters
            subphase_lb, subphase_ub = config.d2o_lb, config.d2o_ub

            # define boolean for whether parameter is fixed or varying
            subphase_vary = True if config.d2o_vary == True else False

        else:
            print(f'Fatal Error: subphase input not acmw or d2o\nsubphase = {self.subphase}')
            sys.exit()


        if self.subphase == 'D2O': math_subphase = 'D_2O'
        if self.subphase == 'ACMW': math_subphase = self.subphase

        str_ = r'$\rho_{'+ math_subphase +'} \ (10^{-6} \AA^2)$'

        # define parameter
        self.rho_subphase = Parameter(subphase_init, bounds=(subphase_lb,subphase_ub), vary=subphase_vary, name=str_)


        # define SLD slab
        self.SLD_subphase = SLD(self.rho_subphase, name=self.subphase)


        # define subphase layer
        if config.rough_bkg_incl == True:
            self.subphase_layer = self.SLD_subphase(0,self.roughness_backing)

        else:
            self.subphase_layer = self.SLD_subphase(0,self.roughness)



    def calc_lipid_SL_SLD(self):

        # caclulate SL of the lipid system
        self.calculatedSL_lip = calc_sld("SL",self.lipid,self.ratio)

        # define SL parameter of the h_lipid tails
        self.b_tails = Parameter(self.calculatedSL_lip.get('tails'), name='SL1')

        # define SL parameter for the head group
        self.b_heads = Parameter(self.calculatedSL_lip.get("head"), name='SL2')


        # calculate lipid SLD
        self.calculatedSLD_lip = calc_sld("SLD",self.lipid,self.ratio)

        # define tails SLD parameter
        # if modelling exchange, then the second contrast is left as a variable
        if config.SLD_exchange_mode == True and '2' in self.contrast:
            self.rho_tails = Parameter(self.calculatedSLD_lip.get('tails'), bounds=(config.sld1_lb,config.sld1_ub), vary=True, name='SLD1 fitted')
            self.rho_tails_2_value = self.rho_tails

            self.rho_heads = Parameter(self.calculatedSLD_lip.get('head'), bounds=(config.sld1_lb,config.sld1_ub), vary=True, name='SLD2 fitted')
            self.rho_heads_2_value = self.rho_heads

        else:
            self.rho_tails = Parameter(self.calculatedSLD_lip.get('tails'), name='SLD1')
            self.rho_tails_1_value = self.rho_tails

            # define heads SLD parameter
            self.rho_heads = Parameter(self.calculatedSLD_lip.get('head'), name='SLD2')
            self.rho_heads_1_value = self.rho_heads



    def define_drug_layer(self):

        # define SLD parameters for hydrogenous and deuterated polyA components
        if 'HPOLYA' in self.contrast:
            self.rho_drug = Parameter(config.rho_h_drug, name='SLD3-h')

        elif 'DPOLYA' in self.contrast:
            self.rho_drug = Parameter(config.rho_d_drug, name='SLD3-d')

        else:
            print('\nFatal Error: Atempted to create drug layer but no polyA in contrast.')
            print(f'Problem contrast = {self.contrast}')
            sys.exit()

        # if not fitting drug layer thickness, constrain it
        if config.with_drug_layer2 == True and config.d3_vary == False:
            self.thickness_layer_3.constraint = config.d3_0 - self.thickness_heads

        # calculate the volume fraction of the drug layer
        self.layer3_vf_drug  = 1 - self.layer3_solvent

        # set as parameter for later constraint
        self.vf_drug_layer3 = Parameter(self.layer3_vf_drug)





    def calc_head_solvent(self):

        # initialise head solvent parameter
        self.head_solvent = Parameter(0.0, name='Solv2')

        # define head volume fraction; headVolFrac = rho1*d1*SL2/rho2*d2*SL1
        self.vf_heads = (self.rho_tails*self.thickness_tails*self.b_heads)/(self.rho_heads*self.thickness_heads*self.b_tails)


        # default head solvent calculation is without drug in head group
        if config.with_drug_layer2 == False:

            # head_solvent = 1 - vf_heads
            self.head_solvent.constraint = 1 - self.vf_heads

        else:

            # head_solvent = 1 - vf_heads - vf_drug_layer3
            self.head_solvent.constraint = 1-self.vf_heads - self.vf_drug_layer3

            # update head SLD; initialise new parameter to prevent circular dependency
            self.rho_heads_merged = Parameter(0, name='SLD2')

            # constrain parameter by averaging head and drug SLDs
            self.rho_heads_merged.constraint = (self.vf_heads*self.rho_heads+self.vf_drug_layer3*self.rho_drug)/(self.rho_heads+self.rho_drug)



    def make_component_slabs(self):

        # create SLD tail object (layer 1)
        self.tails = SLD(self.rho_tails, name='tails')

        # create tails slab
        self.tails_layer = self.tails(self.thickness_tails,self.roughness,self.tail_solvent)


        # create SLD head objects (layer 2)
        if config.with_drug_layer2 == True:
            self.heads = SLD(self.rho_heads_merged, name='heads')

        else:
            self.heads = SLD(self.rho_heads, name='heads')


        # create slab head objects
        self.heads_layer = self.heads(self.thickness_heads,self.roughness,self.head_solvent)



        # create first drug layer (layer 3)
        if config.with_drug_layer3 == True and config.link_drug_thick == False:

            # create SLD object
            self.layer_3 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.third_layer = self.layer_3(self.thickness_layer_3,self.roughness,self.layer3_solvent)

        elif config.with_drug_layer3 == True and config.link_drug_thick == True:

            # create SLD object
            self.layer_3 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.third_layer = self.layer_3(self.linked_drug_thick,self.roughness,self.layer3_solvent)



        # create a second drug layer (layer 4)
        if config.with_drug_layer4 == True and config.link_drug_thick == False:

            # create SLD object
            self.layer_4 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.fourth_layer = self.layer_4(self.thickness_layer_4,self.roughness,self.layer4_solvent)

        elif config.with_drug_layer4 == True and config.link_drug_thick == True:

            # create SLD object
            self.layer_4 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.fourth_layer = self.layer_4(self.linked_drug_thick,self.roughness,self.layer4_solvent)



        # create a third drug layer (layer 5)
        if config.with_drug_layer5 == True and config.link_drug_thick == False:

            # create SLD object
            self.layer_5 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.fifth_layer = self.layer_5(self.thickness_layer_5,self.roughness,self.layer5_solvent)

        elif config.with_drug_layer5 == True and config.link_drug_thick == True:

            # create SLD object
            self.layer_5 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.fifth_layer = self.layer_5(self.linked_drug_thick,self.roughness,self.layer5_solvent)


        # create a fourth drug layer (layer 6)
        if config.with_drug_layer6 == True and config.link_drug_thick == False:

            # create SLD object
            self.layer_6 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.sixth_layer = self.layer_6(self.thickness_layer_6,self.roughness,self.layer6_solvent)

        elif config.with_drug_layer6 == True and config.link_drug_thick == True:

            # create SLD object
            self.layer_6 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.sixth_layer = self.layer_6(self.linked_drug_thick,self.roughness,self.layer6_solvent)


        # create a fifth drug layer (layer 7)
        if config.with_drug_layer7 == True and config.link_drug_thick == False:

            # create SLD object
            self.layer_7 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.seventh_layer = self.layer_7(self.thickness_layer_7,self.roughness,self.layer7_solvent)

        elif config.with_drug_layer7 == True and config.link_drug_thick == True:

            # create SLD object
            self.layer_7 = SLD(self.rho_drug, name='SLD_drug')

            # create drug slab
            self.seventh_layer = self.layer_7(self.linked_drug_thick,self.roughness,self.layer7_solvent)



        # create single 'slab' layer
        if config.surface_excess_mode == True:

            # create thickness parameter for the layer
            self.thickness_slab = Parameter(config.d_slab_0, bounds=(config.d_slab_lb,config.d_slab_ub), vary=config.d_slab_vary, name='Slab Layer Thickness')

            # create SLD object
            self.slab = SLD(self.rho_slab, name='SLD_average_layer')

            # create slab
            self.slab_layer = self.slab(self.thickness_slab,self.roughness,self.slab_solvent)



    def build_structure(self):

        # get experimental data
        exp_data = ReflectDataset(self.file)
        self.exp_list.append(exp_data)

        # check for config mistakes
        if config.with_drug_layer3 == True and config.surface_excess_mode == True:
            print("\nFatal Error: Both drug_layer3 and surface_excess modes selected.\n")
            sys.exit()


        # build structure
        if config.with_drug_layer3 == False and config.with_drug_layer4 == False and config.with_drug_layer5 == False and config.with_drug_layer6 == False and config.with_drug_layer7 == False and config.surface_excess_mode == False:
            structure = self.air | self.tails_layer | self.heads_layer | self.subphase_layer

        elif config.with_drug_layer3 == True and config.with_drug_layer4 == False and config.with_drug_layer5 == False and config.with_drug_layer6 == False and config.with_drug_layer7 == False and config.surface_excess_mode == False:
            structure = self.air | self.tails_layer | self.heads_layer | self.third_layer | self.subphase_layer

        elif config.with_drug_layer3 == True and config.with_drug_layer4 == True and config.with_drug_layer5 == False and config.with_drug_layer6 == False and config.with_drug_layer7 == False and config.surface_excess_mode == False:
            structure = self.air | self.tails_layer | self.heads_layer | self.third_layer | self.fourth_layer | self.subphase_layer

        elif config.with_drug_layer3 == True and config.with_drug_layer4 == True and config.with_drug_layer5 == True and config.with_drug_layer6 == False and config.with_drug_layer7 == False and config.surface_excess_mode == False:
            structure = self.air | self.tails_layer | self.heads_layer | self.third_layer | self.fourth_layer | self.fifth_layer | self.subphase_layer

        elif config.with_drug_layer3 == True and config.with_drug_layer4 == True and config.with_drug_layer5 == True and config.with_drug_layer6 == True and config.with_drug_layer7 == False and config.surface_excess_mode == False:
            structure = self.air | self.tails_layer | self.heads_layer | self.third_layer | self.fourth_layer | self.fifth_layer | self.sixth_layer | self.subphase_layer

        elif config.with_drug_layer3 == True and config.with_drug_layer4 == True and config.with_drug_layer5 == True and config.with_drug_layer6 == True and config.with_drug_layer7 == True and config.surface_excess_mode == False:
            structure = self.air | self.tails_layer | self.heads_layer | self.third_layer | self.fourth_layer | self.fifth_layer | self.sixth_layer | self.seventh_layer | self.subphase_layer


        elif config.surface_excess_mode == True:
            structure = self.air | self.slab_layer | self.subphase_layer

        else:
            print('\nFatal Error: No in-built structure model fits the config instructions.\n')
            sys.exit()




        # store structure data for accessing SLD
        self.struct_list.append(structure)

        # build model from structure
        model = ReflectModel(structure, dq=config.dq, dq_type=config.dq_type)

        # set background
        if self.subphase == 'D2O':
            model.bkg.setp(config.bkg_d2o_0,vary=config.bkg_d2o_vary,bounds=(config.bkg_d2o_lb, config.bkg_d2o_ub))

        elif self.subphase == 'ACMW':
            model.bkg.setp(config.bkg_acmw_0,vary=config.bkg_acmw_vary,bounds=(config.bkg_acmw_lb, config.bkg_acmw_ub))


        # append to list of models
        self.model_list.append(model)

        # generate objective from model and experimental data
        objective = Objective(model, exp_data, use_weights=config.useWeights,)

        # store objective
        self.obj_list.append(objective)




def calcParNum(contrast_list):

    """
    purpose:
        to calculate the number of varied parameters across all the samples

    input:
        contrast_list = list of sample contrasts (NB: each contrast must be different, e.g. _1,_2...)

    analysis:
        if par_vary == True then add parameter
        if par_vary == False then parameter is fixed and not considered

    output:
         nPars = int value
    """

    # initialise number of parameters
    nPars = 0

    # iterate over all contrasts
    for contrast in contrast_list:
        subphase = contrast.split('_')[0]

        # one parameter per instrument background
        if subphase == 'ACMW' and config.bkg_acmw_vary == True:
            nPars += 1

        if subphase == 'D2O' and config.bkg_d2o_vary == True:
            nPars += 1

        # one parameter per subphase SLD
        if subphase == 'ACMW' and config.acmw_vary == True:
            nPars += 1

        if subphase == 'D2O' and config.d2o_vary == True:
            nPars += 1

    # one parameter per thickness layer
    if config.d1_vary == True:
        nPars += 1

    if config.d2_vary == True:
        nPars += 1

    # whether contrasts are fitted with a drug-third layer
    if config.with_drug_layer3 == True:

        if config.d3_vary == True and config.link_drug_thick == False:
            nPars += 1

        if config.solv3_vary == True:
            nPars += 1


    if config.with_drug_layer4 == True:

        if config.d4_vary == True and config.link_drug_thick == False:
            nPars += 1

        if config.solv4_vary == True:
            nPars += 1


    if config.with_drug_layer5 == True:

        if config.d5_vary == True and config.link_drug_thick == False:
            nPars += 1

        if config.solv5_vary == True:
            nPars += 1


    if config.link_drug_thick == True:
        nPars += 1


    # adds parameter for second monolayer SLD
    if config.SLD_exchange_mode == True:
        nPars += 1

    # for cases where backing/background roughness is fitted
    if config.rough_bkg_incl == True and config.rough_bkg_vary == True:
        nPars += 1

    return nPars



def calc_reduced_chisq(exp_list, global_objective, contrast_list):

    """
    purpose:
        to calculate the reduced chi square value, refnx only produces the
        regular chi-sq which does not consider number of points and parameters

    input:
        exp_list = list of experimental data files
        global_objective = the final parameter solution
        contrast_list = list of sample contrasts (NB: each contrast must be different, e.g. _1,_2...)

    output:
        reduced_chisq = float = global_objective.chisqr()/(nPoints-nPars)
    """

    # calculate number of experimental data points
    nPoints = functools.reduce(lambda count, l: count + len(l), exp_list, 0)

    # calculate the number of parameters
    nPars = calcParNum(contrast_list)

    # calculate the reduced chi square
    # global_objective.chisqr()/(nPoints-nPars)
    reduced_chisq = global_objective.chisqr()/(nPoints-nPars)

    # print output in verbose mode
    if config.verbose == True:
        print('\n\nPoints = %d' %nPoints)
        print('nPars = %d' %nPars)
        print('reduced chisqr = global_objective.chisqr()/(nPoints-nPars)')

    return reduced_chisq



def writeParams(parameter_output_str, reduced_chisq, global_objective):

    """
    purpose:
        to store the final parameter solution as a text file for a record of the
        analysis

    input:
        parameter_output_str = string defining path and file save name
        reduced_chisq = float = global_objective.chisqr()/(nPoints-nPars)
        global_objective = the final parameter solution

    output:
        parameters-title.txt file in ../output
    """

    global_objective_file = open(parameter_output_str, 'w')
    global_objective_file.write('\nReduced chiSq = %f\n' %reduced_chisq)
    global_objective_file.write('\nCaclulated from: global_objective.chisqr()/(nPoints-nPars)\n')
    global_objective_file.write(str(global_objective))
    global_objective_file.close()

    return


def fitModel(title, sampleInfo):

    """
    purpose:
        This function does most of the heavy lifting in the program. It takes
        the processed sample information, builds a refnx slab model for each
        sample, fits them, performs MCMC if required, generates the model
        reflectivity data as well as basic parameter analysis.

    input:
        title = the title of the analysis, used for saving outputs
        sampleInfo = pandas dataframe where each row is one sample

    output:
        sampleInfo = updated dataframe to include the model NR data
        global_objective = the final parameter solution
        reduced_chisq = the reduced chi-squared value of the co-refined data
        APM = area per molecule of the co-refined data
    """

    # initialise class
    M = Model()

    # initialise count variable for number of models to cofit
    count_model = 0

    # iterate across each contrast
    for path, contrast, lipid, ratio, label in zip(sampleInfo['filePath'],\
                                                   sampleInfo['contrast'],\
                                                   sampleInfo['lipid'],\
                                                   sampleInfo['ratio'],\
                                                   sampleInfo['label']):

        # define the colour for SLD plots
        colour = config.col_light[count_model]

        # update model count
        count_model +=1

        # load the information into the class
        M.load_contrast(path, contrast, lipid, ratio, title, label, colour)

        # generate the subphase slab
        M.make_backing_slab()

        # calculate the scattering lengths and SLDs of the lipids
        M.calc_lipid_SL_SLD()

        # adds in third layer for drug BUT much match contrast
        # i.e. polyA in contrast and config.with_drug_layer3 = True
        # otherwise chisq will be wrong
        if config.with_drug_layer3 == True:
            M.define_drug_layer()

        # calculate the amount of solvent in the headgroups
        M.calc_head_solvent()

        # generate slabs for each of the components (drug slab deprecated)
        M.make_component_slabs()

        # generate final structure and corresponding objective
        M.build_structure()

    # access completed list of objectives
    obj_list = M.obj_list


    # print line break for debugging clarity
    if config.verbose == True:
        print('\n-------------------------------------------')


    ## Fitting procedure

    # combine all individual objectives into a GlobalObjective
    global_objective = GlobalObjective(obj_list)

    # curve fitter object performs leaste squares fitting
    fitter = CurveFitter(global_objective)

    #constrain t_thick > h_thick
    class DEC(object):
        def __init__(self, pars, global_objective):
            # store the parameters and objective in this object
            # necessary for pickling in the future
            self.pars = pars
            self.objective = global_objective

        def __call__(self, x):
            # update the varying parameters in the objective first
            self.objective.setp(x)
            return float(self.pars[0] - self.pars[1])

    pars = (M.thickness_tails, M.thickness_heads)
    dec = DEC(pars, global_objective)

    thickness_constraint = NonlinearConstraint(dec, 0, np.inf)

    # do fit
    # result, covar : :class:`scipy.optimize.OptimizeResult`, np.ndarray
    # `result.x` contains the best fit parameters
    # `result.covar` is the covariance matrix for the fit.
    # `result.stderr` is the uncertainties on each of the fit parameters.
    fitter.fit(config.algorithm, constraints=(thickness_constraint));

    # calculate reduced chi square value
    reduced_chisq = calc_reduced_chisq(M.exp_list, global_objective, sampleInfo['contrast'])


    ## Markov Chain Monte Carlo
    res = None
    if config.doMCMC == True:

        # MCMC - burn first 400 as initial chain might not be representative of equilibrated system
        fitter.sample(steps=config.MCMC_initSteps, random_state=1, pool=0) # pool = -1 for all cores
        fitter.sampler.reset()

        # MCMC - production run - save 1 in 100 samples to remove autocorrelation, save 15 stes giving 15*200 samples (200 walkers by default)
        res = fitter.sample(config.MCMC_nSteps, nthin=config.MCMC_nThin, random_state=1, pool=0)

        # plot corner plot to show covariance between parameters
        fig = global_objective.corner(show_titles=True)

        # update axis properties
        for ax in fig.get_axes():

            # ax.xaxis.label.set_size(14)
            # ax.yaxis.label.set_size(14)

            # from matplotlib.ticker import FormatStrFormatter
            # ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

            # ax.xaxis.label.set_visible(False)
            # ax.yaxis.label.set_visible(False)


            # ax.xaxis.label.set_weight("bold")
            # ax.yaxis.label.set_weight("bold")

            ax.tick_params(axis='both',
                           direction='in',
                           # labelsize=12,
                           pad=2,
                           bottom=True,
                           top=True,
                           left=True,
                           right=True,
                           )

        plt.savefig('../output/corner-' +title+'.png', # png / svg
            dpi=300,
            bbox_inches='tight')



    # print output if running in verbose mode
    if config.verbose == True:
        print('\n\n\n', global_objective)
        reduced_chisq = calc_reduced_chisq(M.exp_list, global_objective, sampleInfo['contrast'])



    ## Generate model data

    modelQ_list, modelNR_list = [], []

    for idx, model_data in enumerate(M.model_list):

        struct = model_data.structure

        label = sampleInfo['label'][idx]

        col = config.col_light[idx]

        # plot sld
        if config.plotSLD == True:

            # define rc parameters
            mpl.rcParams['axes.linewidth'] = 2
            mpl.rcParams['axes.labelsize'] = 24
            mpl.rcParams['axes.labelsize'] = 24
            mpl.rcParams['xtick.direction'] = 'in'
            mpl.rcParams['ytick.direction'] = 'in'
            mpl.rcParams['xtick.major.width'] = 2
            mpl.rcParams['ytick.major.width'] = 2
            mpl.rcParams['xtick.top'] = True
            mpl.rcParams['ytick.right'] = True
            mpl.rcParams['xtick.labelsize'] = 12
            mpl.rcParams['ytick.labelsize'] = 12
            mpl.rcParams['text.usetex'] = True
            mpl.rcParams['font.family'] = 'serif'
            mpl.rcParams['legend.loc'] = 'best'

            plt.axhline(y=0,color='black',linestyle='-',zorder=0)

            plt.plot(*struct.sld_profile(),lw=3,color=col,label=label)
            plt.xlabel(r'$z$ $\rm{(\mathring{A})}$')
            plt.ylabel(r'$\rho$ $(10^{-6} \rm{\mathring{A}}^{-2})$')
            plt.legend(frameon=False, fontsize=12, ncol=1)

            plt.xlim(-15,100)
            plt.ylim(-0.5,6.5) # 0.25, 2 / -0.5, 6.5

            plt.savefig('../output/sld-'+title+'.png',
                dpi=300,
                bbox_inches='tight')



        # get min and max q values
        qmin, qmax = sampleInfo['qmin'][idx], sampleInfo['qmax'][idx]

        # generate model q data
        q = pd.Series(np.linspace(qmin, qmax, 1001))

        # generate model data
        modelQ_list.append(q)
        modelNR_list.append(pd.Series(model_data(q)))



    # generate h-lipid ACMW contrast
    if config.gen_hlip_acmw_model == True:

        # define model slabs
        sld_air = SLD(0);
        sld_tails = SLD(-0.07299); slab_tails = sld_tails(config.d1_0,config.rough_0,0.0) # h: -0.07299 d: 6.19
        sld_heads = SLD(0.726); slab_heads = sld_heads(config.d2_0,config.rough_0,0.42)
        sld_drug1 = SLD(config.rho_h_drug); slab_drug1 = sld_drug1(config.d3_0,config.rough_0,0.44)
        sld_drug2 = SLD(config.rho_h_drug); slab_drug2 = sld_drug2(config.d4_0,config.rough_0,0.77)
        sld_drug3 = SLD(config.rho_h_drug); slab_drug3 = sld_drug3(config.d5_0,config.rough_0,0.91)
        sld_sub = SLD(0); slab_sub = sld_sub(0,config.rough_0)

        # build structure
        structure = sld_air | slab_tails | slab_heads | slab_drug1 | slab_drug2 | slab_drug3 | slab_sub

        # plot sld
        plt.plot(*structure.sld_profile(),lw=3)
        plt.xlabel(r'$z \rm{[\mathring{A}]}$')
        plt.ylabel(r'$SLD [10^{-6} \rm{\mathring{A}}^{-2}]$')


        # generate model
        model = ReflectModel(structure, dq=config.dq, dq_type=config.dq_type)

        # set background of model
        model.bkg.setp(config.bkg_acmw_0,vary=config.bkg_acmw_vary,bounds=(config.bkg_acmw_lb, config.bkg_acmw_ub))

        # store another q value
        modelQ_list.append(q)
        modelNR_list.append(pd.Series(model_data(q)))

        # this will exceed index of dataframe, therefore get row number and populate
        # current df with None on the new row
        max_df_rows = len(sampleInfo)
        idx = max_df_rows + 1
        sampleInfo.loc[idx] = [None for i in range(len(sampleInfo.loc[0]))]


    # load into dataframe
    sampleInfo['Q_model'] = modelQ_list
    sampleInfo['R_model'] = modelNR_list

    # print(sampleInfo['R_model'])
    # sys.exit()

    ## Calculate APM

    # area per molecule = (SL1/SLD1)/d1 = V1/d1 [Angstroms Sq.]
    # must take SL1 and SLD1 parameters and not volume (V1) directly as this
    # accounts for multicomponent lipid systems
    # does not matter if drug penetrates head layer as calc. only takes tails
    # multiply by 10 for units as in calcPar.py the SLD = SLD * 10 for units
    if config.SLD_exchange_mode == False and config.surface_excess_mode == False:
        APM = 10 * (M.calculatedSL_lip.get('tails')/M.calculatedSLD_lip.get('tails')) / M.thickness_tails.value

    else: APM = 'N/A' # SLD_exchange_mode or surface_excess_mode == True




    ## Calculate charge density

    # taking the last state loaded into M and assuming all contrasts are the same

    # assuming +1e charge per headgroup and APM in angstroms^2
    # if at pH 7.4 then charge number by an additional 0.01 to account for 10 % charge

    # if only one component then force frac=1 to prevent typos in input ratio
    if type(M.ratio) == int:
        mc3_frac = 1

    # if more than one component assume mc3 is written first and convert to frac
    else:
        ratio_list = M.ratio.split(':')
        mc3_frac = float(ratio_list[0]) / 100

    sigma = mc3_frac * (1.609 * 10) / APM



    ## Comparing SLDs

    integral_list = []

    # modelling change in SLD as evidence of lipid exchange
    if config.SLD_exchange_mode == True:

        delta_sld1 = ((M.rho_tails_2_value.value - M.rho_tails_1_value.value) / abs(M.rho_tails_1_value.value) )

        delta_sld2 = ((M.rho_heads_2_value.value - M.rho_heads_1_value.value) / abs(M.rho_heads_1_value.value) )

        print('\nTails SLD fitted for LNP sample.')
        print('Original SLD = %f' %M.rho_tails_1_value.value)
        print('Fitted SLD   = %f' %M.rho_tails_2_value.value)
        print('Percentage change = %f %%' %(delta_sld1*100))

        print('\nHeads SLD fitted for LNP sample.')
        print('Original SLD = %f' %M.rho_heads_1_value.value)
        print('Fitted SLD   = %f' %M.rho_heads_2_value.value)
        print('Percentage change = %f %%' %(delta_sld2*100))

        fig, ax = plt.subplots()

        for idx, struct in enumerate(M.struct_list):

            # access data in array form, x (z pos) then y (SLD)
            struct_data = struct.sld_profile(align=0)

            # assign to x and y variables
            struct_z, struct_SLD = struct_data[0], struct_data[1]
            #print(struct_SLD)

            # quick plot of SLD data
            ax.plot(struct_z, struct_SLD, label=sampleInfo['contrast'][idx])

            # integrate datasets, y then x ordering, should only consider positive values
            I = integrate.simpson(abs(struct_SLD), struct_z)

            # store integral value
            integral_list.append(I)

        # calculate the difference between the two datasets
        delta_integral = ( (integral_list[0] - integral_list[1] ) / (integral_list[0]) )


        print("\nOriginal integrated SLD = %f" %integral_list[0])
        print("New integrated SLD = %f" %integral_list[1])
        print("Percentage change in integrated SLD = %f %%" %(delta_integral*100))

        # set plot parameters
        ax.legend(frameon=False)
        ax.set_ylabel("$\\rho$ / $10^{-6} \AA^{-2}$")
        ax.set_xlabel("z / $\AA$")



    ## Write global objective output to txt file

    if config.writeGlobalObj == True:
        writeParams('../output/par-'+title+'.txt', reduced_chisq, global_objective)


    return sampleInfo, global_objective, reduced_chisq, APM, sigma, res
