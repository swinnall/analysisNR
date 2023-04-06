"""
purpose: 
    user interface file to set parameter values for the fitting procedure   
"""

#########
# Setup #
#########

# constant dQ/Q resolution smearing is employed
dq_type = 'constant'

# dq = 5 gives 5% resolution smearing 
dq = 5

# use experimental uncertainty in calculation of residuals
useWeights = True


#############
# Roughness #
#############

# set roughness for all layers (Alice-EEM: 3.1, Sam-EEM: 3.7, Sam-MC3: 3.5)
rough_0 = 3.5

# set a specific roughness for backing media 
rough_bkg_incl = False 
rough_bkg_0 = 5

# vary the backing media roughness 
rough_bkg_vary = True 
rough_bkg_lb = rough_0
rough_bkg_ub = 30


#########
# Lipid #
#########

# tail thickness (d1) bounds
d1_vary = True  
d1_0  = 14
d1_lb = 10
d1_ub = 30

# head thickness (d2) bounds
d2_vary = False
d2_0  = 5.2
d2_lb = 5
d2_ub = 15

# multiply average chain vol by factor to model lipid chain compaction - CHECK
compact_chains = False
compact_chains_factor = 0.85


########
# Drug #
########

# merge drug with headgroup layer, requires vary_d3 = False; check solv2 is physical
with_drug_layer2 = False

# drug in third layer, must set to False if no drug as effects nPars
with_drug_layer3 = True 

# drug thickness (d3) bounds
d3_vary = True 
d3_0  = 20
d3_lb = 5
d3_ub = 50

# solvent fraction of drug layer
solv3_vary = False
solv3_0  = 0.9
solv3_lb = 0.0
solv3_ub = 1.0


############
# Subphase #
############

# acmw SLD bounds (linked)
acmw_vary = False
acmw_0  = 0.0
acmw_lb = -0.1
acmw_ub = 0.1

# d2o SLD bounds (linked)
d2o_vary = True
d2o_0   = 6.36
d2o_lb  = 6.0
d2o_ub  = 6.36


##############
# Background #
##############

bkg_d2o_vary = False
bkg_d2o_0  = 8e-7
bkg_d2o_lb = 1e-6
bkg_d2o_ub = 1e-5

bkg_acmw_vary = False
bkg_acmw_0  = 8e-7
bkg_acmw_lb = 1e-6
bkg_acmw_ub = 1e-5


############################
# Monte Carlo Markov Chain #
############################

doMCMC = False 

# plot mcmc parameter variation
plotCorner = True if doMCMC == True else False 

# number of burned steps, default = 100 
MCMC_initSteps = 100 # 100 

# number of saved steps, 15*200 samples (200 default walkers)
MCMC_nSteps = 15 # 15  

# how much mcmc chains should be thinned before storing (nThin=1 keeps all values)
MCMC_nThin = 10


#####################
# Exchange Analysis #
#####################

# sets SLD as variable parameter in the second contrast (_2 ending)
# NB: headSolvFrac_2 is calculated from fitted SLD and is not physical
SLD_exchange_mode = False
sld1_lb = 0.1
sld1_ub = 3.0


#####################
# Low Q Composition #
#####################

# fits thickness of a single layer model 
surface_excess_mode = False

d_single_vary = True if surface_excess_mode == True else False
d_single_0 = 22
d_single_lb = 0
d_single_ub = 40

sld_single_0 = 4


#########
# Print #
#########

# levels of printing for debugging
verbose = False
very_verbose = False

# write global objective output to file
writeGlobalObj = True


########
# Plot #
########

# plot sld when making structures
plotSLD = False

# light colours: teal, maroon, orange, olive, purple, silver, yellow, beige
col_light = ['#66B5B3','#C32148','#F38F42','#D4C964',"#CCCCFF",'#CCCCCC','#FFE766','#CCCC99']

# dark colours: teal, maroon, orange, green, purple, charcoal, yellow , indigo
col_dark = ['#116C6E','#8F011B','#C83911','#7D7853',"#3333FF",'#333333','#CCAC00','#273B60']


#########
# Paths #
#########

# user defined sample info database
sample_info_str = 'sample_info.txt'

# input folder - too verbose for samplesNR.txt 
input_sub_list = ['9-13-1016-ILL-MC3','RB2220322-EEM','RB2220338-MC3']
input_sub_str = input_sub_list[0]


#############
# Databases #
#############

# atom coherent scattering lengths [fm], Coh b from https://www.ncnr.nist.gov/resources/n-lengths/
atomSL = {
    "H": -3.739,
    "D": 6.671,
    "C": 6.646,
    "N": 9.36,
    "O": 5.803,
    "P": 5.13,
    "K": 3.67,
    }

# molecular weight lipid database, [g/mol]
lipidMw = {
    "DPPC":            734.039,
    "d-DPPC":          796.43,
    "POPC":            760.07, # 760.076
    "d31-POPC":        791.07, # 791.267
    "POPS":            783.99,
    "Cholesterol":     386.65, # 386.654
    "d45-Cholesterol": 421.8745, # calculated from data in DANIELE/calc_lipids
    "DLin-KC2-DMA":    642.1,
    "DLin-MC3-DMA":    642.09,
    "d62-DLin-MC3-DMA":704.5,
    "Lipid-5":         710.182,
    "d19-Lipid-5":     729.3,
    "DOPE":            744.034,
    "SM":              760.223,
    "LBPA":            792.07,
    "PolyA":           385.31,
    "DMG-PEG-2000":    2509.200,
    }

# chemical structures for each lipid: (struct_head, struct_tail)
lipidStruct = {
    "POPC":            ('N-O8-P-C10-H18', 'C32-H64'), # Yixuan struct. email 03-04-22
    "d31-POPC":        ('N-O8-P-C10-H18', 'C32-D31-H33'),
    "DOPE":            ('N-O8-P-C8-H14', 'C33-H64'),
    "SM":              ('N2-O5-P-C8-H19', 'O1-C33-H64'),
    "LBPA":            ('N-O4-P-C4-H11', 'O6-C38-H71'),
    "Cholesterol":     ('O-H','C27-H45'),
    "d45-Cholesterol": ('O-H','C27-D45'),
    "DLin-MC3-DMA":    ('N-O2-C7-H13', 'C36-H66'),
    "d62-DLin-MC3-DMA":('N-O2-C7-H13', 'C36-H4-D62'),
    "Lipid-5":         ('O-H-N-C7-H14', 'O4-C37-H73'),
    "d19-Lipid-5":     ('O-H-N-C7-H14', 'O4-C37-H54-D19'),
    "DSPC":            ('N-O8-P-C10-H18','C34-H70'),
    "d70-DSPC":        ('N-O8-P-C10-H18','C34-D70'),
    "DMG-PEG-2000":    ('O5-C6-H7','C25-H52'), # polymer: ([O-C2-H4]_44 + O-C3-H7); total: O50-C122-H242
    "PolyA":           ('C10-H13-K-N5-O7-P','H'),
    }

# molecular volumes for each lipid (Angstroms cubed): (head, tail)
lipidMolVol = {
    "POPC":            (344,937),
    "d31-POPC":        (344,937),
    "DOPE":            (236,969),
    "SM":              (274,953),
    "LBPA":            (208,624),
    "Cholesterol":     (5,624),
    "d45-Cholesterol": (5,624),
    "DLin-MC3-DMA":    (260, 1030), # total is 1290 from Arteta paper, 1202 from chemSketch
    "d62-DLin-MC3-DMA":(260, 1030),
    "Lipid-5":         (274, 1000), # total is 1274 from chemSketch
    "d19-Lipid-5":     (274, 1000),
    "DSPC":            (322,1000),
    "d70-DSPC":        (322,1000),
    "DMG-PEG-2000":    (256,767), # From Marianna: DMPE (head 0.25% total vol. 1023) PEG unit = 670
    "PolyA":           (1,1),
    }
