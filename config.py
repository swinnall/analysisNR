"""
purpose:
    user interface file to set parameter values for the fitting procedure

key:
    _vary = set parameter to be fitted (True) or fixed (False)
    _0 = initial value
    _lb = lower boundary
    _ub = upper boundary
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

# algorithm for parameter optimisation (differential_evolution or L-BFGS-B)
algorithm = 'differential_evolution'


#############
# Roughness #
#############

# set roughness for all layers
rough_0 = 3.5

# set a specific roughness for backing media
rough_bkg_incl = False
rough_bkg_0 = 5

# vary the backing media roughness
rough_bkg_vary = False
rough_bkg_lb = rough_0
rough_bkg_ub = 30


#########
# Lipid #
#########

# tails thickness (d1) bounds
d1_vary = True
d1_0  = 15.5
d1_lb = 10
d1_ub = 30

# head thickness (d2) bounds
d2_vary = True
d2_0  = 8
d2_lb = 5
d2_ub = 15

# multiply average chain vol by factor to model lipid chain compaction
compact_chains = False
compact_chains_factor = 0.85


########
# Drug #
########

# SLD of drug (polyA) 10^-6 * A^-2
rho_h_drug = 3.58 # 3.67
rho_d_drug = 4.57 # 4.46


# merge drug with headgroup layer
with_drug_layer2 = False


# force d3 = d4 = d5 (as many as are included)
link_drug_thick = False

d3_link_vary = True 
d3_link_0 = 20
d3_link_lb = 5
d3_link_ub = 50


# drug in third layer
with_drug_layer3 = True

# drug thickness (d3) bounds
d3_vary = True
d3_0  = 20
d3_lb = 5
d3_ub = 50

# solvent fraction of drug layer
solv3_vary = True
solv3_0  = 0.9
solv3_lb = 0.0
solv3_ub = 1.0


# drug in fourth layer
with_drug_layer4 = True

# layer 4 thickness bounds
d4_vary = True
d4_0  = 20
d4_lb = 5
d4_ub = 50

# solvent fraction of layer 4
solv4_vary = True
solv4_0  = 0.9
solv4_lb = 0.0
solv4_ub = 1.0


# drug in 5th layer
with_drug_layer5 = True

# layer 5 thickness bounds
d5_vary = True
d5_0  = 20
d5_lb = 5
d5_ub = 50

# solvent fraction of layer 5
solv5_vary = True
solv5_0  = 0.9
solv5_lb = 0.0
solv5_ub = 1.0



############
# Subphase #
############

# acmw SLD bounds
acmw_vary = False
acmw_0  = 0.0
acmw_lb = -0.1
acmw_ub = 0.1

# d2o SLD bounds
d2o_vary = True
d2o_0   = 6.36
d2o_lb  = 6.0
d2o_ub  = 6.36


##############
# Background #
##############

# d2o backgrounds fitted individually or all fixed to d2o_0 value
bkg_d2o_vary = False
bkg_d2o_0  = 5e-6
bkg_d2o_lb = 1e-6
bkg_d2o_ub = 1e-4

# acmw backgrounds fitted individually or all fixed to acmw_0 value
bkg_acmw_vary = False
bkg_acmw_0  = 9e-6
bkg_acmw_lb = 1e-6
bkg_acmw_ub = 1e-4


############################
# Monte Carlo Markov Chain #
############################

doMCMC = False

# number of burned steps, default = 100
MCMC_initSteps = 100

# number of saved steps, 15*200 samples (200 default walkers)
MCMC_nSteps = 15

# how much mcmc chains should be thinned before storing (nThin=1 keeps all values)
MCMC_nThin = 10


#####################
# Low Q Composition #
#####################

# lowQ is proportional to surface excess & independent of structure (rho*d)
surface_excess_mode = False

# fixed slab SLD
sld_slab_0 = 4

# fits thickness of a single slab model
d_slab_vary = True if surface_excess_mode == True else False
d_slab_0 = 22
d_slab_lb = 0
d_slab_ub = 40


#####################
# Exchange Analysis #
#####################

# sets SLD as variable parameter in the second contrast (_2 ending)
# NB: headSolvFrac_2 is calculated from fitted SLD and is not physical
SLD_exchange_mode = False
sld1_lb = 0.1
sld1_ub = 3.0


#########
# Print #
#########

# levels of printing for debugging
verbose = True

# stop after SLD properties calculated 
SLD_calculator_only = True

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

# size of the experimental datapoints
marker_size = 20

# symbol for the experimental data points
marker_symbol = 'circle'

# thickness of the model line
line_width = 5


#############
# Databases #
#############

# atom coherent scattering lengths [fm], https://www.ncnr.nist.gov/resources/n-lengths/
atomSL = dict(
    H = -3.739,
    D = 6.671,
    C = 6.646,
    N = 9.36,
    O = 5.803,
    P = 5.13,
    K = 3.67,
    )


# lipid parameters
lipid_database = dict(

    # CAPS_LIPID_KEY = dict(mw=total molecular weight [g/mol]
    #                       struct=chemical formula of head and tails
    #                       mvol=molecular volume of head and tails [A^3]
    #                       ),

    CHOL = dict(mw=386.65,
                struct=dict(head='O-H',tails='C27-H45'),
                mvol=dict(head=5,tails=624),
                dfrac=None,
                ),

    D45_CHOL = dict(mw=421.87,
                struct=dict(head='O-H',tails='C27-D45'),
                mvol=dict(head=5,tails=624),
                dfrac=0.79,
                ),
    
    # D45_CHOL = dict(mw=421.87,
    #             struct=dict(head='O-H',tails='C27-D36.34-H9.66'), # 
    #             mvol=dict(head=5,tails=624),
    #             ),

    KC2 = dict(mw=642.10,
                struct=dict(head=None,tails=None),
                mvol=dict(head=None,tails=None),
                dfrac=None,
                ),

    MC3 = dict(mw=642.09,
                struct=dict(head='N-O2-C7-H13',tails='C36-H66'),
                mvol=dict(head=260,tails=1030),
                dfrac=None,
                ),

    D62_MC3 = dict(mw=704.50,
                struct=dict(head='N-O2-C7-H13',tails='C36-H4-D62'),
                mvol=dict(head=260,tails=1030),
                dfrac=None,
                ),

    LIPID_5 = dict(mw=710.18,
                struct=dict(head='O-H-N-C7-H14',tails='O4-C37-H73'),
                mvol=dict(head=274,tails=1000),
                dfrac=None,
                ),

    D19_LIPID_5 = dict(mw=729.30,
                struct=dict(head='O-H-N-C7-H14',tails='O4-C37-H54-D19'),
                mvol=dict(head=274,tails=1000),
                dfrac=None,
                ),

    DPPC = dict(mw=734.04,
                struct=dict(head=None,tails=None),
                mvol=dict(head=None,tails=None),
                dfrac=None,
                ),

    D_DPPC = dict(mw=796.43,
                struct=dict(head=None,tails=None),
                mvol=dict(head=None,tails=None),
                dfrac=None,
                ),

    POPC = dict(mw=760.07,
                struct=dict(head='N-O8-P-C10-H18',tails='C32-H64'),
                mvol=dict(head=344,tails=937),
                dfrac=None,
                ),

    D31_POPC = dict(mw=791.07,
                struct=dict(head='N-O8-P-C10-H18',tails='C32-D31-H33'),
                mvol=dict(head=344,tails=937),
                dfrac=None,
                ),

    POPS = dict(mw=783.99,
                struct=dict(head=None,tails=None),
                mvol=dict(head=None,tails=None),
                dfrac=None,
                ),

    DOPE = dict(mw=744.03,
                struct=dict(head='N-O8-P-C8-H14',tails='C33-H64'),
                mvol=dict(head=236,tails=969),
                dfrac=None,
                ),

    DSPC = dict(mw=None,
                struct=dict(head='N-O8-P-C10-H18',tails='C34-H70'),
                mvol=dict(head=322,tails=1000),
                dfrac=None,
                ),

    D70_DSPC = dict(mw=None,
                struct=dict(head='N-O8-P-C10-H18',tails='C34-D70'),
                mvol=dict(head=322,tails=1000),
                dfrac=None,
                ),

    SM = dict(mw=760.22,
                struct=dict(head='N2-O5-P-C8-H19',tails='O1-C33-H64'),
                mvol=dict(head=274,tails=953),
                dfrac=None,
                ),

    LBPA = dict(mw=792.07,
                struct=dict(head='N-O4-P-C4-H11',tails='O6-C38-H71'),
                mvol=dict(head=208,tails=624),
                dfrac=None,
                ),

    DMG_PEG_2000 = dict(mw=2509.20,
                struct=dict(head='O5-C6-H7',tails='C25-H52'),
                mvol=dict(head=256,tails=767),
                dfrac=None,
                ),

    POLYA = dict(mw=385.31,
                struct=dict(head='C10-H13-K-N5-O7-P',tails='H'),
                mvol=dict(head=1,tails=1),
                dfrac=None,
                ),
    )
