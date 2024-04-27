import chemical as chem_mass

class Molecular(object):
    def __init__(self):
        self.name = ""
        self.chemical = ""
        self.mass = 0
    
    def replace_element(self, replacement):
        '''
        replacemant = [('N','15N'), ('C', '13C'), ...]
        '''
        self.chemical, self.mass = chem_mass.replace_element_and_calc_mass(self.chemical, replacement)

class AminoAcid(Molecular):
    def __init__(self):
        self.name = ""
        self.chemical = ""
        self.mass = 0

class Glyco(Molecular):
    def __init__(self):
        self.name = ""
        self.full_name = ""
        self.chemical = ""
        self.mass = 0
        self.markers = []
        
class Modification(Molecular):
    def __init__(self):
        self.name = ""
        self.short_name = ""
        self.mod_site = ""
        self.chemical = ""
        self.mass = 0
        self.modloss = 0
        
        
class BaseConfig:
    def __init__(self):
        self.gly_dict = {}
        self.mod_dict = {}
        self.aa_dict = {}
    
    def read_aa(self, aaini):
        with open(aaini) as f:
            self.aa_dict = {}
            lines = f.readlines()
            for line in lines:
                if not line.startswith('R'): continue
                aa = AminoAcid()
                items = line.strip().split(" ")
                aa.name = items[0][-1]
                aa.mass = float(items[1])
                aa.chemical = items[2]
                aa.mass = chem_mass.calc_formula_mass(aa.chemical)
                self.aa_dict[aa.name] = aa
                
    def replace_aa_element(self, replacement):
        for key, aa in self.aa_dict.items():
            aa.replace_element(replacement)
    
    def read_mod(self, modini):
        with open(modini) as f:
            self.mod_dict = {}
            lines = f.readlines()
            for line in lines:
                mod = Modification()
                line = line.strip()
                if not line or line.startswith('@') or line.startswith('#') or line.startswith('name'):
                    continue
                items = line.strip().split(' ')
                mod.name = items[0][:items[0].find("=")]
                mod.mod_site = mod.name[mod.name.find('[')+1:mod.name.rfind(']')]
                mod.mass = float(items[2])
                mod.short_name = mod.mod_site + '(' + str(round(mod.mass)) + ')'
                mod.chemical = items[-1]
                mod.mass = chem_mass.calc_formula_mass(mod.chemical)
                if items[4] != '0':
                    mod.modloss = float(items[5])
                    
                self.mod_dict[mod.name] = mod
                
    def replace_mod_element(self, replacement):
        for key, mod in self.mod_dict.items():
            mod.replace_element(replacement)
    
    def read_glyco(self, glycoini):
        with open(glycoini) as f:
            self.gly_dict = {}
            lines = f.readlines()
            for line in lines:
                if not line.startswith('G'): continue
                glyco = Glyco()
                items = line.strip().split(' ')
                glyco.full_name = items[0][items[0].find("=")+1:]
                glyco.name = items[1]
                glyco.mass = float(items[2])
                glyco.chemical = items[4]
                glyco.mass = chem_mass.calc_formula_mass(glyco.chemical)
                if len(items) >= 6 and items[5] != "":
                    markers = items[5].strip(",").split(",")
                    glyco.markers = [float(marker)+1.007276 for marker in markers]
                else:
                    glyco.markers = [glyco.mass+1.007276]
                    
                self.gly_dict[glyco.name] = glyco
                
    def replace_glyco_element(self, replacement):
        for key, glyco in self.gly_dict.items():
            glyco.replace_element(replacement)
                
                