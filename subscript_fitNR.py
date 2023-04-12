from refnx.dataset import ReflectDataset
from refnx.analysis import Objective, Parameter
from refnx.reflect import SLD, ReflectModel
from refnx.analysis import CurveFitter, GlobalObjective
from scipy.optimize import NonlinearConstraint
import matplotlib.pyplot as plt
from scipy import integrate 
import numpy as np
import functools 
import sys

from subscript_calcMonoProp import calc_monolayer_prop
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
        self.thickness_tails = Parameter(config.d1_0, bounds=(config.d1_lb,config.d1_ub), vary=config.d1_vary, name='Tails Thick.')
        self.thickness_heads = Parameter(config.d2_0, bounds=(config.d2_lb,config.d2_ub), vary=config.d2_vary, name='Heads Thick.')
        self.thickness_drug  = Parameter(config.d3_0, bounds=(config.d3_lb,config.d3_ub), vary=config.d3_vary, name='Drug Thick.')
       
        # define tail solvent (always 0, no solvent in tails)
        self.tail_solvent = Parameter(0.0, name='Tail Solv.')
        
        # define how much solvent in drug layer (constant between contrasts)
        if config.with_drug_layer2 == True or config.with_drug_layer3 == True:
            self.drug_solvent = Parameter(config.solv3_0, bounds=(config.solv3_lb,config.solv3_ub), vary=config.solv3_vary, name='Drug Solv.')
        
        # SLD constant (based on reasonable guess) and solvent fixed at 0 
        if config.surface_excess_mode == True:
            self.rho_single = SLD(config.sld_single_0, name='SLD_single_layer')
            self.single_solvent = Parameter(0.0, name='Solv_single_layer')
        
        
        # define lists for data
        self.exp_list = []  
        self.obj_list = []
        self.struct_list = [] 
        self.model_list = []
        
        
    
    
    def load_contrast(self, file, contrast, lipid, ratio):
        
        # define file / sample information parameters 
        self.file = file 
        self.contrast = contrast.split('_')
        
        # receiving a specific lipid and ratio object from main() via lipids.get(i) for e.g.
        self.lipid = lipid
        self.ratio = ratio       
        
        # extract system information from contrast
        self.subphase = self.contrast[0].upper()
        
       
        
        
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

        
        # define parameter 
        self.rho_subphase = Parameter(subphase_init, bounds=(subphase_lb,subphase_ub), vary=subphase_vary, name=self.subphase+' SLD')
        
        
        # define SLD slab 
        self.SLD_subphase = SLD(self.rho_subphase, name=self.subphase)
        
        
        # define subphase layer 
        if config.rough_bkg_incl == True:
            self.subphase_layer = self.SLD_subphase(0,self.roughness_backing) 
        
        else: 
            self.subphase_layer = self.SLD_subphase(0,self.roughness) 
      
    
    
    def calc_lipid_SL_SLD(self):

        # caclulate SL of the lipid system  
        self.calculatedSL_lip = calc_monolayer_prop("SL",self.lipid,self.ratio,1,1)
        
        # define SL parameter of the h_lipid tails 
        self.b_tails = Parameter(self.calculatedSL_lip.get('tails'), name='SL1')
        
        # define SL parameter for the head group 
        self.b_heads = Parameter(self.calculatedSL_lip.get("head"), name='SL2')        


        # calculate lipid SLD 
        self.calculatedSLD_lip = calc_monolayer_prop("SLD",self.lipid,self.ratio,1,1)

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
        if 'hPolyA' in self.contrast:
            self.rho_drug = Parameter(3.67, name='SLD3-h')
            
        elif 'dPolyA' in self.contrast:
            self.rho_drug = Parameter(4.46, name='SLD3-d') 
            
        else:
            print('\nFatal Error: Atempted to create drug layer but no polyA in contrast.\nCheck spelling!')
            sys.exit()

        # if not fitting drug layer thickness, constrain it 
        if config.with_drug_layer2 == True and config.d3_vary == False: 
            self.thickness_drug.constraint = config.d3_0 - self.thickness_heads
         
        # calculate the volume fraction of the drug layer 
        self.vfDrug  = 1 - self.drug_solvent
        
        # set as parameter for later constraint
        self.vf_drug = Parameter(self.vfDrug)



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
            
            # head_solvent = 1 - vf_heads - vf_drug
            self.head_solvent.constraint = 1-self.vf_heads - self.vf_drug
            
            # update head SLD; initialise new parameter to prevent circular dependency
            self.rho_heads_merged = Parameter(0, name='SLD2')
            
            # constrain parameter by averaging head and drug SLDs
            self.rho_heads_merged.constraint = (self.vf_heads*self.rho_heads+self.vf_drug*self.rho_drug)/(self.rho_heads+self.rho_drug)
            
    
    
    
    def make_component_slabs(self):  

        # create SLD head objects 
        if config.with_drug_layer2 == True: 
            self.heads = SLD(self.rho_heads_merged, name='heads')
   
        else:
            self.heads = SLD(self.rho_heads, name='heads')
           
            
        # create slab head objects 
        self.heads_layer = self.heads(self.thickness_heads,self.roughness,self.head_solvent)  


        # create SLD tail object  
        self.tails = SLD(self.rho_tails, name='tail')


        # create tails slab
        self.tails_layer = self.tails(self.thickness_tails,self.roughness,self.tail_solvent) 
            
        
        # create drug layer 
        if config.with_drug_layer3 == True:
            
            # create SLD object 
            self.drug = SLD(self.rho_drug, name='SLD_drug')
            
            # create drug slab
            self.drug_layer = self.drug(self.thickness_drug,self.roughness,self.drug_solvent)
                

        # create single 'average' layer 
        if config.surface_excess_mode == True:
            
            # create thickness parameter for the layer 
            self.thickness_single = Parameter(config.d_single_0, bounds=(config.d_single_lb,config.d_single_ub), vary=config.d_single_vary, name='Single Layer Thickness')            

            # create SLD object 
            self.single = SLD(self.rho_single, name='SLD_average_layer')
            
            # create slab
            self.single_layer = self.single(self.thickness_single,self.roughness,self.single_solvent)



    def build_structure(self):
        
        # get experimental data 
        exp_data = ReflectDataset(self.file)
        self.exp_list.append(exp_data)
        
        # check for config mistakes
        if config.with_drug_layer3 == True and config.surface_excess_mode == True:
            print("\nFatal Error: Both drug_layer3 and surface_excess modes selected.\n")
            sys.exit()
        
        
        # build structure 
        if config.with_drug_layer3 == True: 
            structure = self.air | self.tails_layer | self.heads_layer | self.drug_layer | self.subphase_layer     

        elif config.surface_excess_mode == True:
            structure = self.air | self.single_layer | self.subphase_layer 

        else: 
            structure = self.air | self.tails_layer | self.heads_layer | self.subphase_layer    
        
            
        # plot sld 
        if config.plotSLD == True:
            plt.plot(*structure.sld_profile(),lw=3)
            plt.xlabel(r'$z \rm{[\mathring{A}]}$')
            plt.ylabel(r'$SLD [10^{-6} \rm{\mathring{A}}^{-2}]$')
            
        
        # store structure data for accessing SLD 
        self.struct_list.append(structure)

        # build model from structure 
        model = ReflectModel(structure, dq=config.dq, dq_type=config.dq_type)
        
        # set background 
        if self.subphase == 'D2O': 
            model.bkg.setp(config.bkg_d2o_0,vary=config.bkg_d2o_vary, bounds=(config.bkg_d2o_lb, config.bkg_d2o_ub))

        elif self.subphase == 'ACMW': 
            model.bkg.setp(config.bkg_acmw_0,vary=config.bkg_acmw_vary, bounds=(config.bkg_acmw_lb, config.bkg_acmw_ub))

        
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
        if subphase.upper() == 'ACMW' and config.bkg_acmw_vary == True:
            nPars += 1

        if subphase.upper() == 'D2O' and config.bkg_d2o_vary == True:
            nPars += 1
        
        # one parameter per subphase SLD
        if subphase.upper() == 'ACMW' and config.acmw_vary == True:
            nPars += 1
             
        if subphase.upper() == 'D2O' and config.d2o_vary == True:
            nPars += 1
            
    # one parameter per thickness 
    if config.d1_vary == True:
        nPars += 1
    if config.d2_vary == True:
        nPars += 1
     
    # whether contrasts are fitted with a drug-third layer
    if config.with_drug_layer3 == True:

        if config.d3_vary == True:
            nPars += 1          

        if config.solv3_vary == True:
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
        

def fitModel(title, file_paths, contrast_list, qmin_, qmax_, lipids, ratios):

    """
    purpose: 
        This function does most of the heavy lifting in the program. It takes 
        the processed sample information, builds a refnx slab model for each 
        sample, fits them, performs MCMC if required, generates the model 
        reflectivity data as well as basic parameter analysis. 
        
    input: 
        title = the title of the analysis, used for saving outputs 
        file_paths = dict of file locations by contrast key 
        contrast_list = list of sample contrasts (NB: each contrast must be different, e.g. _1,_2...)
        qmin_ = dict of the minimum q value for a given contrast key 
        qmax_ = dict of the maximum q value for a given contrast key 
        lipids = dict of component lipids string for a given contrast key 
        ratios = dict of component ratios string for a given contrast key 

    output: 
        modelQ = dict of the generated Q vales (x-axis), key is index
        modelNR = dict of the generated R vales (y-axis), key is index
        global_objective = the final parameter solution 
        reduced_chisq = the reduced chi-squared value of the co-refined data 
        APM = area per molecule of the co-refined data 
    """    


    # initialise class 
    M = Model()
    
    # iterate across each contrast 
    for contrast in contrast_list:
        
        # access experiment data file, lipid and ratio
        file, lipid, ratio = file_paths.get(contrast), lipids.get(contrast), ratios.get(contrast)
        
        # load the information into the class 
        M.load_contrast(file, contrast, lipid, ratio)
        
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
    
    
    
    ## Fitting procedure
 
    # combine all individual objectives into a GlobalObjective
    global_objective = GlobalObjective(obj_list)
    
    # curve fitter object performs leaste squares fitting
    fitter = CurveFitter(global_objective)
    
    #constrain t_thick > h_thick 
    class DEC(object):
        def __init__(self, pars, objective):
            # we'll store the parameters and objective in this object
            # this will be necessary for pickling in the future
            self.pars = pars
            self.objective = objective
    
        def __call__(self, x):
            # we need to update the varying parameters in the
            # objective first
            self.objective.setp(x)
            return float(self.pars[0] - self.pars[1])
        
    pars = (M.thickness_tails, M.thickness_heads)
    dec = DEC(pars, global_objective)
    
    thickness_constraint = NonlinearConstraint(dec, 0, np.inf)
     

    # do fit 
    fitter.fit('differential_evolution', constraints=(thickness_constraint));
        
    # calculate reduced chi square value 
    reduced_chisq = calc_reduced_chisq(M.exp_list, global_objective, contrast_list)
 
    
 

    ## Markov Chain Monte Carlo 
           
    if config.doMCMC == True:
        
        # MCMC - burn first 400 as initial chain might not be representative of equilibrated system  
        fitter.sample(steps=config.MCMC_initSteps, random_state=1, pool=-1)
        fitter.sampler.reset()
        
        # MCMC - production run - save 1 in 100 samples to remove autocorrelation, save 15 stes giving 15*200 samples (200 walkers by default)
        res = fitter.sample(config.MCMC_nSteps, nthin=config.MCMC_nThin, random_state=1, pool=-1)
        
        # plot corner plot to show covariance between parameters 
        global_objective.corner()
        plt.savefig('../output/corner-' +title+'.png',
            format='png',
            dpi=400,
            bbox_inches='tight')    
        


    # print output if running in verbose mode 
    if config.verbose == True:
        print('\n\n\n', global_objective)
        reduced_chisq = calc_reduced_chisq(M.exp_list, global_objective, contrast_list)



    ## Generate model data 
    
    modelQ  = {}
    modelNR = {}
    
    for idx, model_data in enumerate(M.model_list):
        
        # get relevant contrast 
        contrast = contrast_list[idx]
        
        # get min and max q values 
        qmin, qmax = qmin_.get(contrast), qmax_.get(contrast)
        
        # generate model q data
        q = list(np.linspace(qmin, qmax, 1001)) 
        
        # generate model data 
        modelQ[idx] = q
        modelNR[idx] = list(model_data(q))
    
    
    
    ## Calculate APM
    
    # area per molecule = (SL1/SLD1)/d1 = V1/d1 [Angstroms Sq.]
    # must take SL1 and SLD1 parameters and not volume (V1) directly as this
    # accounts for multicomponent lipid systems 
    # does not matter if drug penetrates head layer as calc. only takes tails
    # multiply by 10 for units as in calcPar.py the SLD = SLD * 10 for units 
    if config.SLD_exchange_mode == False and config.surface_excess_mode == False:
        APM = 10 * (M.calculatedSL_lip.get('tails')/M.calculatedSLD_lip.get('tails')) / M.thickness_tails.value

    else: APM = 'N/A' # SLD_exchange_mode or surface_excess_mode == True

    
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
            ax.plot(struct_z, struct_SLD, label=contrast_list[idx])
   
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
        writeParams('../output/parameters-'+title+'.txt', reduced_chisq, global_objective)
        
        
        
    return modelQ, modelNR, global_objective, reduced_chisq, APM