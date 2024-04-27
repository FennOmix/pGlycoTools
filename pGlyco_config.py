from base_config import BaseConfig
import sys
        
class pGlycoConfig:
    def __init__(self, pGlyco_cfg = None):
        self.glycan_type = "N-Glycan"
        self.fragmentation = "HCD"
        self.glyco_names = "H,N,A,F".split(",")
        self.core_glycans = []
        self.core_formulas = []
        self.glyco_var_mod = []
        self.max_glyco_mod = 0
        
        self.protein_fix_mod = ["Carbamidomethyl[C]"]
        self.protein_var_mod = ["Oxidation[M]"]
        
        self.isppm = 1
        self.tol = 20
        self.pre_isppm = 1
        self.pre_tol = 10
        
        self.self_boosted = True
        self.iteration = 10
        self.output_dir = "."
        self.spec_file_type = "mgf"
        self.mgf_list = []
        self.glycoini = "glyco.ini"
        self.modini = "modification.ini"
        self.aaini = "aa.ini"
        self.fasta = ""
        self.enzyme = ""
        # self.enzyme = "Trypsin KR _ C _ _ N"
        
        self.MS1_RT_win = 120*2
        self.MS1_Scan_win = 400
        
        self.aa_dict = {}
        self.gly_dict = {}
        self.mod_dict = {}
        self.base_conf = BaseConfig()
        
        self.verbose = False
        self.pGlycoType = "pGlycoDB"
        
        self.FDR = 0.01
        self.percolator = 0
        self.FMM_for_peptide_FDR = 0
        if pGlyco_cfg: self.read_cfg(pGlyco_cfg)
        
    def result_template(self):
        return self.pGlycoType + "-GP-Raw%d.txt"
    
    def read_aa(self, aaini):
        self.base_conf.read_aa(aaini)
        self.aa_dict = self.base_conf.aa_dict
    
    def read_mod(self, modini):
        self.base_conf.read_mod(modini)
        self.mod_dict = self.base_conf.mod_dict
    
    def read_glyco(self, glycoini):
        self.base_conf.read_glyco(glycoini)
        self.gly_dict = self.base_conf.gly_dict
        
    def __str__(self):
        s = "pGlyco config:\n"
        s += "  glycan type = %s\n"%self.glycan_type
        s += "  fragmentation = %s\n"%self.fragmentation
        s += "  output_dir = %s\n"%self.output_dir
        s += "  glycoini = %s\n"%self.glycoini
        s += "  modini = %s\n"%self.modini
        return s

    def read_cfg(self, pGlyco_cfg):
        cfg = {}
        with open(pGlyco_cfg) as f:
            lines = f.readlines()
            for line in lines:
                items = line.strip().split("=")
                if len(items) < 2: continue
                cfg[items[0]] = items[1]
            self.mgf_list = []
            for i in range(len(lines)):
                if lines[i].startswith("spectrum_total="):
                    self.mgf_list = []
                    n_spec = int(lines[i].strip().split("=")[1])
                    for j in range(i+1, i+n_spec+1):
                        self.mgf_list.append(lines[j].strip().split("=")[1])
                    break
            for i in range(len(lines)):
                if lines[i].startswith("msms_file_path="):
                    self.mgf_list.append(lines[i].strip().split("=")[1])
                    
        if "spec_file_type" in cfg: self.spec_file_type = cfg['spec_file_type']
        if "msms_file_type" in cfg: self.spec_file_type = cfg['msms_file_type']
        
        if 'fasta' in cfg: self.fasta = cfg['fasta']
        
        if 'enzyme' in cfg: 
            self.enzyme = cfg['enzyme']
            if "digestion" in cfg and cfg['digestion'].startswith('semi'):
                self.enzyme = ""
        
        if "pGlyco_type" in cfg:
            self.pGlycoType = cfg["pGlyco_type"]
            
        if 'aaini' in cfg:
            self.aaini = cfg['aaini']
            self.read_aa(self.aaini)
        else:
            self.aaini = 'aa.ini'
            self.read_aa(self.aaini)
                    
        if "glycoini" in cfg:
            self.glycoini = cfg["glycoini"]
            self.read_glyco(self.glycoini)
        else:
            self.glycoini = "glyco.ini"
            self.read_glyco(self.glycoini)
        
        if "modini" in cfg:
            self.modini = cfg["modini"]
            self.read_mod(self.modini)
        else:
            self.modini = "modification.ini"
            self.read_mod(self.modini)
            
        
        if "fragment_tolerance" in cfg:
            self.tol = float(cfg["fragment_tolerance"])
        
        if "fragment_tolerance_type" in cfg:
            if cfg["fragment_tolerance_type"].lower() == "ppm":
                self.isppm = 1
            else:
                self.isppm = 0
                
        if "precursor_tolerance" in cfg:
            self.pre_tol = float(cfg["precursor_tolerance"])
        
        if "precursor_tolerance_type" in cfg:
            if cfg["precursor_tolerance_type"].lower() == "ppm":
                self.pre_isppm = 1
            else:
                self.pre_isppm = 0
        
        if "glycan_type" not in cfg:
            self.glycan_type = "N-Glycan"
        else:
            self.glycan_type = cfg["glycan_type"]
        
        if "glycan_core" not in cfg:
            if self.glycan_type == "N-Glycan":
                glycan_core = "N(1),N(2),N(2)H(1),N(2)H(2),N(2)H(3),N(1)F(1),N(2)F(1)".split(",")
            else:
                glycan_core = ["N(1)"]
        else:
            glycan_core = cfg["glycan_core"].split(",")
        self.core_formulas = glycan_core
        
        if "glycan_fix_mod" in cfg:
            fixmod = cfg["glycan_fix_mod"].split(",")
            if fixmod[0] != "":
                fixmod = [mod.split("~") for mod in fixmod]
                for m1, m2 in fixmod:
                    m1 += "("
                    m2 += "("
                    for i in range(len(self.core_formulas)):
                        if self.core_formulas[i].startswith(m1):
                            self.core_formulas[i] = self.core_formulas[i].replace(m1, m2)
                
        if "glycan_var_mod" in cfg:
            self.glyco_var_mod = cfg["glycan_var_mod"].split(",")
            self.glyco_var_mod = [mod.split("~") for mod in self.glyco_var_mod]
            
        if "max_var_mod_on_glycan" in cfg:
            self.max_glyco_mod = int(cfg["max_var_mod_on_glycan"])
            
        if "protein_fix_mod" in cfg:
            self.protein_fix_mod = []
            for mod in cfg["protein_fix_mod"].strip(',').split(','):
                if mod: self.protein_fix_mod.append(mod)
        if "protein_var_mod" in cfg:
            self.protein_var_mod = []
            for mod in cfg["protein_var_mod"].strip(',').split(','):
                if mod: self.protein_var_mod.append(mod)
            
        if "fragmentation_type" in cfg:
            self.fragmentation = cfg["fragmentation_type"]
            
        if "output_dir" in cfg:
            self.output_dir = cfg["output_dir"]
        else:
            print("[ERROR] no output_dir in cfg file")
            sys.exit(-1)
        
        if 'percolator' in cfg:
            self.percolator = int(cfg['percolator'])
            
        if 'MS1_RT_win' in cfg:
            self.MS1_RT_win = float(cfg['MS1_RT_win'])
            
        if 'MS1_Scan_win' in cfg:
            self.MS1_RT_win = int(cfg['MS1_Scan_win'])
            
        if 'FDR' in cfg:
            self.FDR = float(cfg['FDR'])
            
        if 'FMM_for_peptide_FDR' in cfg:
            self.FMM_for_peptide_FDR = int(cfg['FMM_for_peptide_FDR'])

    def set_glyco_names(self, glyco_names):
        self.glyco_names = glyco_names
        
if __name__ == "__main__":
    def write_default():
        with open("append_this_to_pGlyco_cfg.txt","w") as f:
            f.write("pGlyco_type=pGlycoDB\n")
            f.write("percolator=1\n")
            f.write("FDR=0.01\n")
            f.write("[Quant]\n")
            f.write("MS1_RT_win=240\n")
            f.write("MS1_Scan_win=400\n")
    write_default()