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

    output: 
        getMonolayerPar method returns SL, SLD and head solvent fraction
    """        

    def __init__(self, lipidNames, molRatios):

        # read system parameters
        self.lipidNames  = lipidNames
        self.molRatios   = molRatios

        # import databases
        self.atomSL        = config.atomSL
        self.lipidDatabase = config.lipid_database 

        # initialise new variables
        self.normMolRatios      = []
        self.lipidSL            = {}
        self.sumAvSL            = 0
        self.sumAvVol           = 0
        self.sumAvSLD           = 0
        self.totalLipidVol      = {'head': 0, 'tails': 0}
        self.avVol              = {'head': 0, 'tails': 0}
        self.avSL               = {'head': 0, 'tails': 0}
        self.avSLD              = {'head': 0, 'tails': 0}


    # converts molar ratio txt input to a normalised version e.g. 4:5 ==> 4/9:5/9
    def normaliseMolarRatios(self):
        
        totMol = 0
        for i, value in enumerate(self.molRatios):
            
            try: totMol += float(value)
            
            except ValueError: 
                print(f'\nFatal Error: cannot convert molar ratio to float\nmolar ratio = {self.molRatios}')
                sys.exit()

        for i, value in enumerate(self.molRatios):
            self.normMolRatios.append( float(value) / totMol )

        if config.verbose == True:
            print("\nLipid names:\n%s\n\nInput molar ratios:\n%s" %(self.lipidNames, self.molRatios))
            print("\nNormalised molar ratios:\n%s" %self.normMolRatios)




    # calculates the scattering length of each lipid component
    def calcSL(self):

        for i, lipid in enumerate(self.lipidNames):
            self.lipidSL[lipid] = {'head': 0, 'tails': 0}

            for j, struct in enumerate(['head','tails']):

                ## SUM THE SCATTERING LENGTHS OF THE ATOMS                 

                # splits head/tails into list of constituent atoms
                splitStruct = re.split('-', self.lipidDatabase.get(lipid).get('struct').get(struct))

                # this loop iterates across ele_strments in a given head/tails and sums the atomic scattering lengths
                for ele_str in splitStruct:
                    
                    # split into ['atom','number of given atom']
                    ele_split_str = list(filter(None, re.split(r'(\d+)', ele_str)))
                    
                    # get atom string 
                    atom_str = ele_split_str[0]
                    
                    # get number of atoms, if only one char in list then atom_num = 1 
                    if len(ele_split_str) == 1:
                        atom_num = 1
                    else:
                        atom_num = int(ele_split_str[1])

    
                    if atom_str == 'D' and self.lipidDatabase.get(lipid).get('dfrac') != None:
                        
                        # get deuteration coeficient from config database 
                        dfrac = self.lipidDatabase.get(lipid).get('dfrac')
                        
                        # calculate the associated scattering length 
                        self.lipidSL[lipid][struct] += dfrac * self.atomSL.get(atom_str) * atom_num
                        
                        # now must account for remaining hydrogenous fraction 
                        hfrac = 1 - dfrac 
                        
                        # change atom string to H and multiply by hfrac
                        self.lipidSL[lipid][struct] += hfrac * self.atomSL.get('H') * atom_num
                        
                        if config.verbose == True:
                            print(f'\nAccounted for {dfrac*100} %% deuteration')
                            print(f'number of D = {dfrac:.3} * {atom_num} = {dfrac*atom_num:.3}')
                            print(f'number of H added = {hfrac:.3} * {atom_num} = {hfrac*atom_num:.3} ')
                    
                    # no deuteration fraction applied 
                    else: 
                        self.lipidSL[lipid][struct] += self.atomSL.get(atom_str) * atom_num                        
                
                
                # multiply total SL of a ipid's head/tails by corresponding normalised molar ratio 
                # convert from fm to m
                self.avSL[struct] += (10**-15) * self.normMolRatios[i] * self.lipidSL[lipid][struct]


        # calculate the average SL of the whole averaged lipid
        self.sumAvSL = self.avSL.get("head") + self.avSL.get("tails")
        
        
        if config.verbose == True:
             print("\nlipidSL [fm]:\n%s" %self.lipidSL)
             print("\navSL:\n%s" %self.avSL)
             print("\nsumAvSL:\n%.3e" %self.sumAvSL)




    def calcVol(self):

        for i, lipid in enumerate(self.lipidNames):

            for j, struct in enumerate(['head','tails']):
                
                # calculate total lipid volumes # convert from A^3 to m
                self.avVol[struct] += (10**-30) * self.normMolRatios[i] * self.lipidDatabase.get(lipid).get('mvol').get(struct)   


        # after average volumes calculated, account for chain compaction if needed
        if config.compact_chains == True:
            self.avVol['tails'] = config.compact_chains_factor * self.avVol['tails']
        
        
        # calculate the average volume of the whole averaged lipid
        self.sumAvVol = self.avVol.get("head") + self.avVol.get("tails")
        

        if config.verbose == True:
            print("\navVol:\n%s" %self.avVol)
            print("\nsumAvVol:\n%.3e" %self.sumAvVol)



    # calculates the scattering length density of each lipid component
    def calcSLD(self):
        
        # head and tail SLDs 
        self.avSLD['tails'] = (10**-14)*self.avSL['tails']/self.avVol['tails']
        self.avSLD['head'] = (10**-14)*self.avSL['head']/self.avVol['head']

        # calculate the average SLD of the whole molecule
        # note that: SL_t/V_t + SL_h/V_h != (SL_t + SL_h)/(V_t + V_h)
        # therefore calculate from (SL_t + SL_h)/(V_t + V_h)
        self.sumAvSLD = (10**-14)*(self.sumAvSL / self.sumAvVol)
     
        if config.verbose == True:
            print(f'\navSLD["tails"] = avSL_tails/avVolv_tails = {self.avSLD["tails"]:.3f} 10^-6 A^-2')
            print(f'avSLD["head"] = avSL_head/avVol_head = {self.avSLD["head"]:.3f} 10^-6 A^-2')
            print(f'\ntotal_avSLD = sumAvSL/sumAvVol = {self.sumAvSLD:.3f} 10^-6 A^-2')

  



# function to check whether a string contains a number
def hasNumbers(inputString):
    return bool(re.search(r'\d', inputString))




def calc_sld(calcType,membrane,lipidRatio):

    """
    purpose: 
        function that builds the Monolayer class for calculating scattering and
        head solvent fraction values 
        
    input: 
        calcType = the type of analysis intended (SL, SLD or head_vol)
        membrane = the string of lipids for a given sample 
        lipidRatio = the string of ratios for each lipid component 
    
    analysis:
        Membrane class 
    
    output: 
        scattering length, scattering length density
    """    
    
    # print line break for debugging clarity 
    if config.verbose == True:
        print('\n-------------------------------------------')
        print('calcType = %s' %calcType)
    
    # get list of lipids within membrane with ratio
    lipids = re.split(':',membrane)
    ratios = re.split(':',str(lipidRatio))

    # create class instance with input variables
    m = Monolayer(lipids, ratios)

    # converts 3:5 to 3/8:5/8
    m.normaliseMolarRatios()

    # calculate coherent scattering lengths
    m.calcSL()

    
    if calcType.upper() == 'SL':
        
        # update SL to units of fm 
        m.avSL.update((x, y*(10**15)) for x, y in m.avSL.items())
        
        return m.avSL

    if calcType.upper() == 'SLD':
        
        # calculate average lipid structure volumes
        m.calcVol()

        # divide scattering length by the molecular volume
        m.calcSLD()
        
        # stop after first SLD calculation to check values
        if config.SLD_calculator_only == True:
            sys.exit()
        
        return m.avSLD

