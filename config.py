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

# name of the file containing sample information
sample_info_fname = 'sample_info_chol.txt'

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
rough_0 = 3.1

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
d1_0  = 20
d1_lb = 5
d1_ub = 30

# head thickness (d2) bounds
d2_vary = True
d2_0  = 5.76
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
d3_vary = False
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
d4_vary = False
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
d5_vary = False
d5_0  = 20
d5_lb = 5
d5_ub = 50

# solvent fraction of layer 5
solv5_vary = True
solv5_0  = 0.9
solv5_lb = 0.0
solv5_ub = 1.0


# drug in 6th layer
with_drug_layer6 = False

# layer 6 thickness bounds
d6_vary = False
d6_0  = 20
d6_lb = 5
d6_ub = 50

# solvent fraction of layer 6
solv6_vary = True
solv6_0  = 0.9
solv6_lb = 0.0
solv6_ub = 1.0


# drug in 7th layer
with_drug_layer7 = False

# layer 7 thickness bounds
d7_vary = False
d7_0  = 20
d7_lb = 5
d7_ub = 50

# solvent fraction of layer 7
solv7_vary = True
solv7_0  = 0.9
solv7_lb = 0.0
solv7_ub = 1.0


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
bkg_d2o_vary = True
bkg_d2o_0  = 1.2e-5 # 1e-6
bkg_d2o_lb = 1e-7
bkg_d2o_ub = 1e-4

# acmw backgrounds fitted individually or all fixed to acmw_0 value
bkg_acmw_vary = True
bkg_acmw_0  = 1.2e-5 # 1e-6
bkg_acmw_lb = 1e-7
bkg_acmw_ub = 1e-4


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

# levels of printing for debugging or SL calculations
verbose = False

# stop after SLD properties calculated
SLD_calculator_only = False

# write global objective output to file
writeGlobalObj = True


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


########
# Plot #
########

# plot sld when making structures
plotSLD = True

# plot only the experimental data (no fitting of the model)
plotOnlyExp = False

# generate and plot an h-lipid acmw model dataset
gen_hlip_acmw_model = False

# axis ranges
xrange = [0, 0.35] # low q: [0.02,0.04], full q: [0, 0.35]
yrange = [-6.2,0.2] # acmw: [-7.2,-3.2], d2o: [-7.2,0.2], log range by exponent: 10^0=1, 10^5=100000

# light colours: blue, green, red, brown
col_light = ['#1e91ff','#6EC531','#fe0000','#fbc48d','#fc8eac']

# dark colours: blue, green, red, brown
col_dark = ['#00008a','#00693C','#8b0000','#8a4704','#c154c1']

# size of the experimental datapoints
marker_size = 10

# marker opacity
opacity = 1.0

# symbol for the experimental data points
#symbol_list = ['circle-open', 'square-open', 'triangle-up-open', 'diamond-open']
symbol_list = ['circle', 'square', 'triangle-up', 'diamond', 'circle-open']
marker_symbol = symbol_list

# thickness of the model line
line_width = 3

# define the title of the legend
legend_title = 'Apr. 2023: pH 6.5' # 'MC3:chol (56.5:43.5), pH 7.4<br>0.08 mg/ml PolyA'

# define legend position
legend_x = 0.50 # prev: 0.35
legend_y = 1.0

# define legend font size
legend_fs = 22



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
    #                       dfrac=percentage deuteration as a fraction, None==1 is default
    #                       ),

# EEM LIPIDS

    POPC = dict(mw=760.07,
                struct=dict(head='N-O8-P-C11-H20',tails='C31-H62'),
                mvol=dict(head=319,tails=944),
                dfrac=None,
                # ref - both volumes: 10.1016/j.bpj.2014.09.027
                # ref - total mvol=1242: 10.1021/la501835p
                ),

    D31_POPC = dict(mw=791.07,
                struct=dict(head='N-O8-P-C11-H20',tails='C31-D31-H31'),
                mvol=dict(head=319,tails=944),
                dfrac=None,
                ),

    DOPE = dict(mw=744.03,
                struct=dict(head='N-O8-P-C8-H14',tails='C33-H64'),
                mvol=dict(head=312,tails=923), # values 0.25 % of total volume?
                dfrac=None,
                # ref - total vol 1235: 10.1016/S0006-3495(97)78067-6
                ),

    SM = dict(mw=760.22,
                struct=dict(head='N2-O6-P-C9-H19',tails='C32-H64'),
                mvol=dict(head=274,tails=937),
                dfrac=None,
                # ref - both vols but not sure which tails 10.1021/acs.jpcb.0c03389
                ),

    CHOL = dict(mw=386.65,
                struct=dict(head='O-H',tails='C27-H45'),
                mvol=dict(head=29,tails=601),
                dfrac=None,
                # ref - head vol by reported mass density (table 1): 10.1016/j.bpj.2013.07.003
                # ref - chol at 630 +/- 10: 10.1016/j.chemphyslip.2006.04.002
                # ref - popc and chol volumes gives mvol chol 655 A^3: 10.1021/la501835p
                ),

    D45_CHOL = dict(mw=421.87,
                struct=dict(head='O-H',tails='C27-D45'),
                mvol=dict(head=29,tails=601),
                dfrac=0.79,
                ),

# OTHER MEMB. LIPIDS

    POPS = dict(mw=783.99,
                struct=dict(head=None,tails=None),
                mvol=dict(head=None,tails=None),
                dfrac=None,
                ),

    LBPA = dict(mw=792.07,
                struct=dict(head='N-O4-P-C4-H11',tails='O6-C38-H71'),
                mvol=dict(head=208,tails=624),
                dfrac=None,
                ),

# LNP LIPIDS

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

    DMG_PEG_2000 = dict(mw=2509.20,
                struct=dict(head='O5-C6-H7',tails='C25-H52'),
                mvol=dict(head=256,tails=767),
                dfrac=None,
                ),

# CATIONIC IONISABLE LIPIDS

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
                dfrac=0.996842,
                ),

# MISC

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

    POLYA = dict(mw=385.31,
                struct=dict(head='C10-H13-K-N5-O7-P',tails='H'),
                mvol=dict(head=1,tails=1),
                dfrac=None,
                ),

    NA = dict(mw=1,
                struct=dict(head='H',tails='H'),
                mvol=dict(head=1,tails=1),
                dfrac=None,
                ),
    )
