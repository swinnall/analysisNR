import re
import sys
import config

class Monolayer:

    """
    purpose: 
        Monolayer is a class that builds on itself to calculate the average 
        properties of the components within the same 
        
    input: 
        lipidNames = list of lipid component names 
        molRatios = list of molar % amounts of the respective lipids 
        thickness = dictionary of head and tail thicknesses
        monolayerPar = tuple of dicts (needed if adding another component to the monolayer)

    output: 
        getMonolayerPar method returns SL, SLD and head solvent fraction
        
    comment:
        this class is overqualified for the fitting procedure and could be 
        reduced in complexity in the future 
    """        

    def __init__(self, lipidNames, molRatios, thickness, monolayerPar):

        # read system parameters
        self.lipidNames  = lipidNames
        self.molRatios   = molRatios
        self.thickness   = thickness
        self.nLipids     = len(self.lipidNames)

        # unpack monolayer parameters (0s if 1st run; each tuple ele is dict[struct])
        self.monolayerMolVol = monolayerPar[0]
        self.monolayerSL     = monolayerPar[1]
        self.monolayerSLD    = monolayerPar[2]

        # import databases
        self.lipidStruct = config.lipidStruct
        self.atomSL      = config.atomSL
        self.lipidMolVol = config.lipidMolVol

        # initialise new variables
        self.headVolFrac        = 0
        self.twoSolv            = 0
        self.d3                 = 0
        self.threeSolv          = 0
        self.twoSLD_H2O         = 0
        self.twoSLD_D2O         = 0
        self.drugOverlapSLD_H2O = 0
        self.drugOverlapSLD_D2O = 0
        self.normMolRatios      = []
        self.lipidSL            = {}
        self.sumAvSL            = {}
        self.sumAvSLD           = {}
        self.totalLipidVol      = {'head': 0, 'tails': 0}
        self.avLipidVol         = {'head': 0, 'tails': 0}
        self.avSL               = {'head': 0, 'tails': 0}
        self.avSLD              = {'head': 0, 'tails': 0}


    # converts molar ratio txt input to a normalised version e.g. 4:5 ==> 4/9:5/9
    def normaliseMolarRatios(self):
        
        totMol = 0
        for i, value in enumerate(self.molRatios):
            totMol += float(value)

        for i, value in enumerate(self.molRatios):
            self.normMolRatios.append( float(value) / totMol  )

        if config.verbose == True and "Monolayer" not in self.lipidNames:
            print("\nLipid names:\n%s\n\nInput molar ratios:\n%s" %(self.lipidNames, self.molRatios))

        if config.very_verbose == True:
            print("\nNormalised molar ratios:\n%s" %self.normMolRatios)


    # calculates the total lipid volume
    def calcTotalLipidVol(self):

        for i, lipid in enumerate(self.lipidNames):

            # check lipid exists in database
            if self.lipidStruct.get(lipid) == None and lipid == "Monolayer": pass
            elif self.lipidStruct.get(lipid) == None:
                print("\nError: Lipid type not found in Lipid Molecular Formula Database.")
                print("Lipid: %s" %lipid)
                sys.exit()

            # calculate total lipid volumes
            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.monolayerMolVol[struct]
                else:
                    self.totalLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerMolVol['head']  = self.totalLipidVol.get('head')
            self.monolayerMolVol['tails'] = self.totalLipidVol.get('tails')

        if config.very_verbose == True:
            print('\nTotal Volume:\n%s' %self.totalLipidVol)




    # calculates the scattering length of each lipid component
    def calcSL(self):

        for i, lipid in enumerate(self.lipidNames):
            self.lipidSL[lipid]       = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.lipidSL[lipid][struct] = self.monolayerSL.get(struct)

                else:
                    # splits head/tail into list of constituent atoms
                    splitStruct = re.split('-', self.lipidStruct.get(lipid)[j])

                    # this loop iterates across elements in a given head/tail and sums the atomic scattering lengths
                    for ele in splitStruct:

                        # multiply scattering length of atom by number of atoms
                        if hasNumbers(ele) == True:

                            # x is split into ['atom','number of given atom']
                            x = list(filter(None, re.split(r'(\d+)', ele)))

                            self.lipidSL[lipid][struct] += self.atomSL.get(x[0]) * int(x[1])


                        # add the scattering length of the single atom identified
                        elif hasNumbers(ele) == False:
                            self.lipidSL[lipid][struct] += self.atomSL.get(ele[0])


                # multiply total lipid scattering length of a given lipid's head/tail by corresponding vol frac
                self.avSL[struct] += self.normMolRatios[i] * self.lipidSL[lipid][struct]


        # calculate the average SL of the whole averaged lipid
        self.sumAvSL = self.avSL.get("head") + self.avSL.get("tails")

        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerSL['head']  = self.avSL.get('head')
            self.monolayerSL['tails'] = self.avSL.get('tails')

        if config.very_verbose == True:
            print("\nLipid scattering lengths:\n%s" %self.lipidSL)

        if config.verbose == True:
            print("\nAverage SL:\n%s" %self.avSL)

        if config.very_verbose == True:
            print("\nSummed average SL:\n%s" %self.sumAvSL)


    def calcAvLipidVol(self):

        for i, lipid in enumerate(self.lipidNames):

            # calculate total lipid volumes
            for j, struct in enumerate(['head','tails']):

                if lipid == "Monolayer":
                    self.avLipidVol[struct] += self.normMolRatios[i] * self.monolayerMolVol[struct]
                else:
                    self.avLipidVol[struct] += self.normMolRatios[i] * self.lipidMolVol.get(lipid)[j]


        # if accounting for chain compaction
        if config.compact_chains == True:
            compact_chains_factor = config.compact_chains_factor
            self.avLipidVol['tails'] = compact_chains_factor * self.avLipidVol['tails']


        if config.verbose == True:
            print("\nAverage Lipid Head/Tail Volume:\n%s" %self.avLipidVol)



    # calculates the scattering length density of each lipid component
    def calcSLD(self):

        # take avSL and divide by avStructVol (tails, heads)
        for i, struct in enumerate(['head','tails']):
            self.avSLD[struct] = 10 * self.avSL[struct] / self.avLipidVol[struct]

        # calculate the average SLD of the whole averaged lipid
        self.sumAvSLD = self.avSLD.get("head") + self.avSLD.get("tails")

        # save Monolayer monolayer struct volumes on the first iteration
        if "Monolayer" not in self.lipidNames:
            self.monolayerSLD['head']  = self.avSLD.get('head')
            self.monolayerSLD['tails'] = self.avSLD.get('tails')

        if config.verbose == True:
            print("\nAverage SLD:\n%s" %self.avSLD)

        if config.very_verbose == True:
            print("\nSummed average SLD:\n%s" %self.sumAvSLD)


    # calculates volume fraction of the head group based on SLD
    def calcHeadVolumeFraction(self):

        # call SL
        SL1 = self.avSL.get("tails")
        SL2 = self.avSL.get("head")

        # call SLD
        rho1 = self.avSLD.get("tails")
        rho2 = self.avSLD.get("head")

        # call membrane thickness
        if "Monolayer" not in self.lipidNames:
            d1 = self.thickness.get("tails")
            d2 = self.thickness.get("head")

        # set new membrane thickness if intended
        elif "Monolayer" in self.lipidNames:
            if config.updateMonolayerThickness == True:
                d1 = config.new_d1
                d2 = config.new_d2
                print("\nYou have selected to update the monolayer thickness such that d1 = %.3f and d2 = %.3f." %(d1,d2))

            else:
                d1 = self.thickness.get("tails")
                d2 = self.thickness.get("head")
                print("\nMonolayer thicknesses are as in instructions file; d1 = %.3f and d2 = %.3f." %(d1,d2))

        # calculate volume fraction
        self.headVolFrac = (rho1 * d1 * SL2 ) / ( rho2 * d2 * SL1)
        self.twoSolv     = (1-self.headVolFrac)*100

        # solvent volume in head group
        if config.verbose == True:
            print("\nHead volume fraction: %f" %self.headVolFrac)
            print("\n2-solv = %f" %self.twoSolv)

        return (1-self.headVolFrac)

    # method that returns monolayer information for adding extra lipid components
    def getMonolayerPar(self):
        return (self.monolayerMolVol, self.monolayerSL, self.monolayerSLD)



# function to check whether a string contains a number
def hasNumbers(inputString):
    return bool(re.search(r'\d', inputString))




def calc_monolayer_prop(calcType,membrane,lipidRatio,t_thick,h_thick):

    """
    purpose: 
        function that builds the Monolayer class for calculating scattering and
        head solvent fraction values 
        
    input: 
        calcType = the type of analysis intended (SL, SLD or head_vol)
        membrane = the string of lipids for a given sample 
        lipidRatio = the string of ratios for each lipid component 
        t_thick = tail thickness (only needed for head_solv calculation)
        h_thick = head thickness (only needed for head_solv calculation)
    
    analysis:
        Membrane class 
    
    output: 
        scattering length, scattering length density or head solvent fraction
    """    


    # get list of lipids within membrane with ratio
    lipids = re.split(':',membrane)
    ratios = re.split(':',str(lipidRatio))

    # structure thickness information into dict
    thickness = {
        'head':  float(h_thick),
        'tails': float(t_thick)
    }

    # initialise monolayer information
    monolayerMolVol = {'head': 0, 'tails': 0}
    monolayerSL     = {'head': 0, 'tails': 0}
    monolayerSLD    = {'head': 0, 'tails': 0}
    monolayerPar    = (monolayerMolVol, monolayerSL, monolayerSLD)


    # create class instance with input variables
    m = Monolayer(lipids, ratios, thickness, monolayerPar)

    # converts 3:5 to 3/8:5/8
    m.normaliseMolarRatios()

    # calculates the total lipid volume of the monolayer (needed for addLipidToMonolayer)
    m.calcTotalLipidVol()

    # calculate coherent scattering lengths
    m.calcSL()

    # calculate average lipid structure volumes
    m.calcAvLipidVol()

    # divide scattering length by the molecular volume
    m.calcSLD()

    # calculate volume fraction of the headgroups
    solvFrac = m.calcHeadVolumeFraction()
    
    # get monolayer parameters: (self.monolayerMolVol, self.monolayerSL, self.monolayerSLD)
    monolayerPar = m.getMonolayerPar()
    
    if calcType.upper() == 'SL':
        return monolayerPar[1]

    if calcType.upper() == 'SLD':
        return monolayerPar[2]
    
    if calcType.upper() == 'solvFrac':
        return solvFrac
