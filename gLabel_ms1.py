#coding=utf-8
'''
Created on 2013.12.13

@author: dell
'''

import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import wx
import copy
import datetime, time, sys
import traceback
from scipy.stats import pearsonr

# from ms_reader import GetMSReader
import AAMass
from isotope import IsotopeSimple as isotp

import spectrum_index as spi
from mass_utils import GlyIonCalc, PepIonCalc
from base_config import BaseConfig
from pGlyco_config import pGlycoConfig


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

isotope_len = 11

# labels = "|N~15N|C~13C"
labels = ""
# if len(sys.argv) > 1:
#     labels = sys.argv[1]
# if labels.lower() == "none":
#     labels = ""
# elif labels == "15N":
#     labels = "|N~15N"
# elif labels == "13C":
#     labels = "|C~13C"
# elif labels == "13C15N" or labels == "15N13C":
#     labels = "|N~15N|C~13C"
    
label_names = labels.split("|")
for i in range(len(label_names)):
    if not label_names[i]:
        label_names[i] = "Unlabel"
labels = [[tuple(one_label.split("~")) for one_label in label.split(";")] if label else [] for label in labels.split("|")]

colors = ["red","blue","orange"]

base_molecular_list = []
isotope_dist_list = []
for i in range(len(labels)):
    base_molecular = BaseConfig()
    base_molecular.read_aa('aa.ini')
    base_molecular.read_mod('modification.ini')
    base_molecular.read_glyco('glyco.ini')
    base_molecular.replace_aa_element(labels[i])
    base_molecular.replace_glyco_element(labels[i])
    base_molecular_list.append(base_molecular)
    
    conf = pGlycoConfig()
    conf.mod_dict = base_molecular.mod_dict
    conf.aa_dict = base_molecular.aa_dict
    conf.gly_dict = base_molecular.gly_dict
    isotope_dist_list.append(isotp(conf, isotope_len))

pepcalc_list = [PepIonCalc(base_molecular_list[i].mod_dict, base_molecular_list[i].aa_dict) for i in range(len(labels))]
glycalc_list = [GlyIonCalc(['H','N','A','G','F'], base_molecular_list[i].gly_dict) for i in range(len(labels))]

def reset_glycalc_list(glyco_names):
    for i in range(len(labels)):
        glycalc_list[i] = GlyIonCalc(glyco_names, base_molecular_list[i].gly_dict)
        
def GetModList(modstr):
    if modstr == "" or modstr.lower() == "null": return []
    modinfo = modstr.strip(';').split(";")
    modlist = []
    for mod in modinfo:
        site, modname = mod.split(",")
        site = int(site)
        modlist.append((site, modname))
    modlist.sort(key = lambda x : x[0])
    return modlist

class Marker:
    def __init__(self, glyco = "H", shape = "o", color = "green", marker = "", alt = "white", fill = "full"):
        self.glyco = glyco
        self.shape = shape
        self.color = color
        self.alt = alt
        self.fill = fill
        if marker == "": self.oxonium_ions = []
        else: self.oxonium_ions = [eval(i) for i in marker.strip().split(";")]

fillstyles = ["right","top","left","bottom"]
def fill_rotation(fill, rotation = 90):
    rotation = int(rotation/90)
    if fill == "full" or fill == "none": return fill
    elif fill == "right": return fillstyles[rotation%4]
    elif fill == "top": return fillstyles[(1+rotation)%4]
    elif fill == "left": return fillstyles[(2+rotation)%4]
    elif fill == "bottom": return fillstyles[(3+rotation)%4]
    else: return fill

class CGlobalConfig:
    def __init__(self):
        self.glycandb = "./pGlyco.gdb"
        self.glycan_type = 'N'
        self.activation_type = "HCD"
        self.ResultHasFragments = 1
        self.plotDecoyPeptide = 1
        self.plotDecoyGlycan = 0
        self.plotMaxGlycanFDR = 1
        self.plotMaxPeptideFDR = 1
        self.plotMaxTotalFDR = 1
        self.plotMinGlycanScore = 0
        self.plotMinPeptideScore = 0
        self.plotMinTotalScore = 0
        self.isPlot = 1
        self.isBatchPlot = 0
        self.glyco_as_text = 0
        self.max_oxonium_mz = 400
        self.glycoini = ""
        self.modini = ""
        self.dpi = 120
        self.SetDefault()
        self.ReadConfig()
        
    def SetDefault(self):
        self.used_markers = []
        self.used_markers.append(Marker('H', 'o', 'green', '145.0495347452;163.0600994315;366.1394719645'))
        self.used_markers.append(Marker('N', 's', 'blue', '138.0552587690;168.0655191604;186.0760838467;204.0866485330'))
        self.used_markers.append(Marker('A', 'D', 'purple', '274.0921278414;292.1026925277;657.2348884922'))
        self.used_markers.append(Marker('G', 'D', 'cyan', '290.08704246349997;308.0976071497999;673.2298031143'))
        self.used_markers.append(Marker('F', '^', 'red', '147.0651848094;350.1445573424'))
        self.used_glyco = ['H','N', 'A', 'G', 'F']
        # reset_glycalc_list(self.used_glyco)
        # self.used_shape = ["o","s","D","D","^"]
        # self.used_color = ["g","b","purple","cyan","r"]
        # self.used_alt = ["white"]*5
        # self.used_fill = ["full"]*5
        # self.used_marker = [[145.0495347452, 163.0600994315, 366.1394719645],    #Hex
                              # [138.0552587690, 168.0655191604, 186.0760838467, 204.0866485330],    #HexNAc
                              # [274.0921278414, 292.1026925277, 657.2348884922],    #NeuAc
                              # [290.08704246349997, 308.0976071497999, 673.2298031143],    #NeuGc
                              # [146.0579088094+aamass.mass_proton, 146.0579088094+204.0866485330]]    #dHex
        self.GetCoreFrags()
        
    def num_format(self, glycocomp):
        return "(" + ",".join([str(g) for g in glycocomp]) + ")"
        
    def short_format(self, glycocomp):
        if sum(glycocomp) == 0: return ""
        ret = ""
        for i in range(len(glycocomp)):
            if glycocomp[i] > 0: ret += "%s(%d)"%(self.used_glyco[i],glycocomp[i])
        return ret
    
    def FormatGlycan(self, glycocomp):
        return self.short_format(glycocomp)
        
    def GetCoreFrags(self):
        if self.glycan_type == 'N': self.GetNLinkedCoreFrags()
        else: self.GetOtherCoreFrags()
        
    def GetOtherCoreFrags(self):
        self.core_frags = []
        self.core_fuc = []
    
    def ExtendCoreFrags(self, glycan, gly_list):
        for g in self.core_frags:
            if np.all(np.array(g) <= np.array(glycan)):
                gly_list.append(g)
        try:
            idx_F = self.used_glyco.index('F')
            if glycan[idx_F] > 0: 
                for g in self.core_fuc:
                    if np.all(np.array(g) <= np.array(glycan)):
                        gly_list.append(g)
        except:
            pass
    
    def GetNLinkedCoreFrags(self):
        self.core_frags = []
        self.core_fuc = []
        try:
            idx_N = self.used_glyco.index('N')
            Y1 = [0]*len(self.used_glyco)
            Y1[idx_N] = 1
            self.core_frags.append(tuple(Y1))
            Y2 = [0]*len(self.used_glyco)
            Y2[idx_N] = 2
            self.core_frags.append(tuple(Y2))
            try:
                idx_H = self.used_glyco.index('H')
                Y3 = copy.deepcopy(Y2)
                Y3[idx_H] = 1
                self.core_frags.append(tuple(Y3))
                Y4 = copy.deepcopy(Y2)
                Y4[idx_H] = 2
                self.core_frags.append(tuple(Y4))
                Y5 = copy.deepcopy(Y2)
                Y5[idx_H] = 3
                self.core_frags.append(tuple(Y5))
            except:
                pass
            try:
                idx_F = self.used_glyco.index('F')
                Y2_fuc = copy.deepcopy(Y1)
                Y2_fuc[idx_F] = 1
                self.core_fuc.append(tuple(Y2_fuc))
                Y3_fuc = copy.deepcopy(Y2)
                Y3_fuc[idx_F] = 1
                self.core_fuc.append(tuple(Y3_fuc))
            except:
                pass
        except:
            pass
    
    def UseGlyco(self, glyco):
        old_glyco = copy.deepcopy(config.used_glyco)
        self.RearrangeGlyco(glyco)
        return old_glyco
    
    def RearrangeGlyco(self, glyco):
        new_glyco = glyco
        new_markers = []
        # new_col = []
        # new_marker = []
        # new_alt = []
        # new_fill = []
        for g in glyco:
            try:
                idx = self.glyco.index(g)
                new_markers.append(self.marker_list[idx])
                # new_col.append(self.glyco_col[idx])
                # new_marker.append(self.oxonium_markers[idx])
            except:
                print('{} is not a glyco unit in glabel.gconf (units = {})'.format(g, ','.join(self.glyco)))
                print('please add {} into glabel.gconf, and then click "Load gLabel config"'.format(g))
                
        self.used_glyco = new_glyco
        self.used_markers = new_markers
        reset_glycalc_list(self.used_glyco)
        # self.used_shape = new_shape
        # self.used_color = new_col
        # self.used_marker = new_marker
        self.GetCoreFrags()
        
    def ReadConfig(self, conf_file = "glabel.gconf"):
        f = open(conf_file)
        lines = f.readlines()
        f.close()
        
        self.marker_list = []
        self.glyco = []
        # self.glyco_shape = []
        # self.glyco_col = []
        # self.oxonium_markers = []
        
        def glyco_shape_color(line):
            items = line.split("=")
            g = items[0].strip()
            items = items[1].split(",")
            items = [item.split(":") for item in items]
            kargs = dict([(key.strip(),val.strip()) for key, val in items])
            kargs['glyco'] = g
            self.marker_list.append(Marker(**kargs))
                
            # items = items[1][items[1].find("shape:")+len("shape:"):items[1].find(",color")].strip()
            # color = items[1][items[1].find("color:")+len("color:"):items[1].find(",marker")].strip()
            # markers = items[1][items[1].find("marker:")+len("marker:"):].strip()
            # markers = [float(marker) for marker in markers.split(',') if marker != ""]
            
            self.glyco.append(g)
            # self.glyco_shape.append(shape)
            # self.glyco_col.append(color)
            # self.oxonium_markers.append(markers)
            
        
        for line in lines:
            if line.startswith("#"): continue
            elif line.startswith("glycandb"):
                self.glycandb = line[line.find("=")+1:].strip()
            elif line.startswith("glycan_type"):
                self.glycan_type = line[line.find("=")+1:].strip()
            elif line.startswith("result_has_fragments"):
                self.ResultHasFragments = int(line[line.find("=")+1:].strip())
            elif line.startswith("plot_decoy_peptide"):
                self.plotDecoyPeptide = int(line[line.find("=")+1:].strip())
            elif line.startswith("plot_decoy_glycan"):
                self.plotDecoyGlycan = int(line[line.find("=")+1:].strip())
            elif line.startswith("plot_max_glycan_FDR"):
                self.plotMaxGlycanFDR = float(line[line.find("=")+1:].strip())
            elif line.startswith("plot_max_peptide_FDR"):
                self.plotMaxPeptideFDR = float(line[line.find("=")+1:].strip())
            elif line.startswith("plot_max_total_FDR"):
                self.plotMaxTotalFDR = float(line[line.find("=")+1:].strip())
            elif line.startswith("plot_min_glycan_score"):
                self.plotMinGlycanScore = float(line[line.find("=")+1:].strip())
            elif line.startswith("plot_min_peptide_score"):
                self.plotMinPeptideScore = float(line[line.find("=")+1:].strip())
            elif line.startswith("plot_min_total_score"):
                self.plotMinTotalScore = float(line[line.find("=")+1:].strip())
            elif line.startswith("is_batch_plot"):
                self.isBatchPlot = int(line[line.find("=")+1:].strip())
            elif line.startswith("glyco_as_text"):
                self.glyco_as_text = int(line[line.find("=")+1:].strip())
            elif line.startswith("activation_type"):
                self.activation_type = line[line.find("=")+1:].strip()
            elif line.startswith("save_dpi"):
                self.dpi = int(line[line.find("=")+1:].strip())
            elif line.startswith("glycoini"):
                self.glycoini = line[line.find("=")+1:].strip()
            elif line.startswith("modini"):
                self.modini = line[line.find("=")+1:].strip()
            elif line.startswith("AA"):
                items = line[line.find("=")+1:].strip().split(":")
                aamass.aa_mass_dict[items[0].strip()] = float(items[1].strip())
            elif "shape:" in line and "color:" in line:
                glyco_shape_color(line)
        self.RearrangeGlyco(self.used_glyco)
        # if self.glycoini: aamass.__read_glyco_ini__(self.glycoini)
        # if self.modini: aamass.__read_mod_ini__(self.modini)

config = CGlobalConfig()
  
fontsize = 12
markersize = 8

vfactor = 0.025

#end parameters

def CalcIsotopeDist(gpsm, isotope_dist):
    return isotope_dist.get_distribution(gpsm.peptide, gpsm.mod, gpsm.glycan, config.used_glyco)
#
def CalcGlycoPeptideMz(gpsm, pepcalc, glycalc):
    modmass = pepcalc.calc_mod_mass_list(gpsm.peptide, gpsm.mod)
    bions = pepcalc.calc_b_ions(gpsm.peptide, modmass)
    pepmass = pepcalc.calc_pepmass_from_b(gpsm.peptide, modmass, bions)
    glymass = glycalc.calc_glycan_mass(gpsm.glycan)
    return (pepmass + glymass)/gpsm.charge + pepcalc_list[0].base_mass.mass_proton
    
class gPSM:
    def __init__(self, spec="", peptide="", glycan=None, glycan_list=[]):
        self.spec = spec
        self.peptide = peptide
        self.glycan = glycan
        self.charge = 0
        self.glycan_list = glycan_list
        self.glycan_decoy = False
        self.mod = ""
        self.modlist = []
        self.glysite = -1
        self.scan = 0

class gPSMList:
    def __init__(self, db_file = config.glycandb):
        self.GlycanCol = 'Glycan'
        self.psmlist = {}
        
    def ReadDenovoRes(self, psm_file):
        self.psmlist = {}
        f = open(psm_file)
        line = f.readline()
        items = line[:-1].split("\t")
        item_idx = {}
        for i in range(len(items)):
            item_idx[items[i]] = i
            if items[i].startswith('Glycan('):
                self.GlycanCol = items[i]
                glystr = self.GlycanCol[self.GlycanCol.find('(')+1:self.GlycanCol.find(')')]
                config.RearrangeGlyco(glystr.split(','))
        for i in range(len(items)):
            item_idx[items[i]] = i
        while True:
            line = f.readline()
            if line == "": break
            items = line.split("\t")
            if items[item_idx["Rank"]] != "1": continue
            if config.plotDecoyPeptide == 0 and items[item_idx["GlyDecoy"]] == "1": continue
            if config.plotDecoyPeptide == 0 and items[item_idx["PepDecoy"]] == "1": continue
            if "GlycanFDR" in item_idx and config.plotMaxGlycanFDR != 1 and float(items[item_idx["GlycanFDR"]]) > config.plotMaxGlycanFDR: continue
            if "PeptideFDR" in item_idx and config.plotMaxPeptideFDR != 1 and float(items[item_idx["PeptideFDR"]]) > config.plotMaxPeptideFDR: continue
            if "GlyScore" in item_idx and float(items[item_idx["GlyScore"]]) < config.plotMinGlycanScore: continue
            if "PepScore" in item_idx and float(items[item_idx["PepScore"]]) < config.plotMinPeptideScore: continue
            if "TotalScore" in item_idx and float(items[item_idx["TotalScore"]]) < config.plotMinTotalScore: continue
            gpsm = gPSM()
            gpsm.spec = items[item_idx["PepSpec"]]
            glyspec = items[item_idx["GlySpec"]]
            gpsm.scan = int(items[item_idx['Scan']])
            glycos = items[item_idx[self.GlycanCol]].strip().split(" ")
            gpsm.glycan = tuple(int(glyco) for glyco in glycos)
            
            gpsm.peptide = items[item_idx["Peptide"]]
            gpsm.glysite = int(items[item_idx["GlySite"]]) - 1
            if "Charge" in item_idx: gpsm.charge = int(items[item_idx["Charge"]])
            else: gpsm.charge = int(glyspec.split('.')[-3])
            gpsm.mod = items[item_idx["Mod"]]
            
            self.psmlist[gpsm.spec] = gpsm
            if gpsm.spec != glyspec:
                gpsm = copy.deepcopy(gpsm)
                gpsm.spec = glyspec
                self.psmlist[gpsm.spec] = gpsm
        f.close()

#
class Label(object):
    '''
    classdocs
    '''
#     def __del__(self):
#         if self.output_info:
#             self.outmsg.close()
    def __init__(self, tol, tol_type="Da"):
        '''
        Constructor
        '''
        self.tol = tol
        self.tol_type = tol_type
        self.reader = None
        self.gpsms = gPSMList()
        self.output_info = False
        self.plot_peptide = True
        self.plot_glycan = True
        self.show_mass = False
        self.RT_win = 120
        
        self.max_plot_mz = 2100.0
        
    def ReadRAW(self, input_file):
        self.input_spec = os.path.split(input_file)[1]
        if self.reader is not None: self.reader.close()
        pf1 = input_file[:-3]+'pf1'
        if os.path.isfile(input_file):
            self.reader = GetMSReader(input_file)
        elif os.path.isfile(pf1):
            self.reader = GetMSReader(pf1)
            
    def FindMS1PeaksByMS2Scan(self, scan):
        scan -= 1
        while scan > 0:
            if scan in self.reader.scanidx:
                RT = self.reader.scanidx[scan][1]
                peaklist = self.reader.read_a_peaklist(scan)
                return peaklist, RT, scan
            scan -= 1
        return None, None, None
    
    def SeeOnePlot_new(self, gpsm):
        mz_ints, RT, ms1_scan = self.FindMS1PeaksByMS2Scan(gpsm.scan)
        if mz_ints is not None:
            print("%s: %s-%s, MS1Scan=%d, RT=%.3f"%(gpsm.spec, gpsm.peptide, config.FormatGlycan(gpsm.glycan), ms1_scan, RT))
        else:
            print("cannot find ms1 scan for '%s'" %gpsm.spec)
            return False
        
        
        ############### glycan and peptide mass ###############
        ions = np.array([])
        ion_types = []
        ion_colors = np.array([])
        iso_dist = []
        for i in range(len(labels)):
            mz = CalcGlycoPeptideMz(gpsm, pepcalc_list[i], glycalc_list[i])
            dist,mono = CalcIsotopeDist(gpsm, isotope_dist_list[i])
            label_ion_types = ["" for i in range(isotope_len)]
            label_ion_types[mono] = label_names[i]
            label_ions = mz*np.ones(isotope_len)+pepcalc_list[i].base_mass.mass_isotope/gpsm.charge*(np.arange(isotope_len)-mono)
            ions = np.append(ions, label_ions)
            ion_types.extend(label_ion_types)
            ion_colors = np.append(ion_colors, np.ones(isotope_len, dtype=int)*i)
            iso_dist.append(dist)
        ion_colors = np.int32(ion_colors)
        ion_types = np.array(ion_types)
        # print(ions)
        ############### end glycan and peptide mass ###############
        
        min_x = np.min(ions)-10
        max_x = np.max(ions)+20
            
        ############### init plot ###############
        xmz = mz_ints[:,0]
        yint = mz_ints[:,1]
        
        yint = yint[np.logical_and(xmz <= max_x, xmz >= min_x)]
        xmz = xmz[np.logical_and(xmz <= max_x, xmz >= min_x)]
        
        max_inten = np.max(yint)
        
        mz_ints = np.append(xmz.reshape(-1,1), yint.reshape(-1,1), axis=1)
        
        if config.isPlot:
            gs = gridspec.GridSpec(4, 1, height_ratios=[1.5,5,1,3.5])
            self.fig = plt.figure(figsize=(16,10)) #
            ax3 = self.fig.add_subplot(gs[0,0])
            ax1 = self.fig.add_subplot(gs[1,0])
            ax2 = self.fig.add_subplot(gs[2,0])
            axLC = self.fig.add_subplot(gs[3,0])
            max_height = max_inten * 1.6
            ax1.vlines([0,max_x],[0.2,0],[max_height, max_height],color="w")
            ax1.hlines(0,0,max_x,linewidth=0.5)
        
            #ax1.set_ylim([0,max_inten])
            ax1.vlines(xmz, [0]*len(yint), yint, color="gray")
        ############### end init plot ###############
        
        def plot_LC(axLC, ms1_scan, ms1_RT):
            ms1_scan_list = [(ms1_scan, ms1_RT)]
            scan = ms1_scan - 1
            while scan > 0:
                if scan in self.reader.scanidx:
                    RT = self.reader.scanidx[scan][1]
                    if abs(RT-ms1_RT) > self.RT_win: break
                    ms1_scan_list.append((scan, RT))
                scan -= 1
            scan = ms1_scan + 1
            while scan <= self.reader.last_scan:
                if scan in self.reader.scanidx:
                    RT = self.reader.scanidx[scan][1]
                    if abs(RT-ms1_RT) > self.RT_win: break
                    ms1_scan_list.append((scan, RT))
                scan += 1
            ms1_scan_list.sort()
                
            ions = np.array([])
            ion_types = []
            ion_colors = np.array([])
            iso_dist = []
            for i in range(len(labels)):
                mz = CalcGlycoPeptideMz(gpsm, pepcalc_list[i], glycalc_list[i])
                dist,mono = CalcIsotopeDist(gpsm, isotope_dist_list[i])
                label_ion_types = ["" for i in range(isotope_len)]
                label_ion_types[mono] = label_names[i]
                label_ions = mz*np.ones(isotope_len)+pepcalc_list[i].base_mass.mass_isotope/gpsm.charge*(np.arange(isotope_len)-mono)
                ions = np.append(ions, label_ions)
                ion_types.extend(label_ion_types)
                ion_colors = np.append(ion_colors, np.ones(isotope_len, dtype=int)*i)
                iso_dist.append(dist)
            ion_colors = np.int32(ion_colors)
            ion_types = np.array(ion_types)
            
            extracted_ions = [[] for i in range(len(ions))]
            LC_max_inten = 0
            for _scan, _RT in ms1_scan_list:
                peaklist = self.reader.read_a_peaklist(_scan)
                idx, mass_tol, intens = self.MatchPeak(peaklist, ions)
                for i in range(len(ions)):
                    if idx[i] != 1:
                        extracted_ions[i].append((_RT, intens[i], colors[ion_colors[i]]))
                _max_inten = max(intens)
                if LC_max_inten < _max_inten: LC_max_inten = _max_inten
            axLC.plot([ms1_RT, ms1_RT], [0, LC_max_inten], "--", color = "gray")
                        
            for i in range(len(extracted_ions)):
                if len(extracted_ions[i]) > 0:
                    # print(extracted_ions[i])
                    axLC.plot([items[0] for items in extracted_ions[i]], [items[1] for items in extracted_ions[i]],
                        "-" if ion_types[i] else "--", color = extracted_ions[i][0][-1], linewidth = (2 if ion_types[i] else 0.5))
                
                
        idx, mass_tol, intens = self.MatchPeak(mz_ints, ions)
        idx = np.array(idx)
        mass_tol = np.array(mass_tol)
        
        intens = np.array(intens)
        R_dict = {}
        for i in range(len(labels)):
            R = pearsonr(iso_dist[i], intens[(i*isotope_len):((i+1)*isotope_len)])[0]
            R_dict[label_names[i]] = R
            print("{}: Pearson R = {}".format(label_names[i], R))
        
        ions = ions[idx != -1]
        ion_types = ion_types[idx != -1]
        ion_colors = ion_colors[idx != -1]
        mass_tol = mass_tol[idx != -1]
        idx = idx[idx != -1]

        if config.isPlot:
        
            plot_LC(axLC, ms1_scan, RT)
            
            ########## plot
            order = np.argsort(idx)
            mass_tol = mass_tol[order]
            ions = ions[order]
            ion_types = ion_types[order]
            ion_colors = ion_colors[order]
            idx = idx[order]
            peptide = gpsm.peptide
            plotpeptide = peptide
            peplen = len(peptide)
            plotmod = "Mod: "
            
            modlist = []
            if gpsm.mod:
                mods = gpsm.mod.strip(";").split(";")
                for mod in mods:
                    site, mod = mod.split(",")
                    modlist.append((int(site), mod))
                modlist.sort()
                for modidx, modname in modlist:
                    if modidx == 0: plotmod += "NTerm+%d"%round(base_molecular_list[0].mod_dict[modname].mass) + ";"
                    elif modidx == peplen+1: plotmod += "CTerm+%d"%round(base_molecular_list[0].mod_dict[modname].mass) + ";"
                    else: plotmod += gpsm.peptide[modidx-1] + str(modidx) + ('+' if base_molecular_list[0].mod_dict[modname].mass > 0 else '') + str(round(base_molecular_list[0].mod_dict[modname].mass)) + ";"
            else:
                plotmod = "noMod"
            
            xmztol = xmz
            xmz = xmz[idx]
            yint = yint[idx]
            
            text_offset = 0
            fontsize = 8
            markersize = 8

            height = 0.02
            baseheight = max_inten+height*max_inten
                    
            for i in range(len(ions)):
                if ion_types[i]:
                    ax1.plot([xmz[i],xmz[i]],[yint[i], baseheight], "--", color="gray", linewidth=0.5)
                    ax1.text(xmz[i],baseheight, "%s: R=%.2f\nm/z=%.2f"%(ion_types[i],R_dict[ion_types[i]],xmz[i]),
                        rotation=90, fontsize=fontsize, horizontalalignment="center",verticalalignment="bottom")
                ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=colors[ion_colors[i]], linewidth=2)
                ax2.plot(xmz[i], mass_tol[i], color=colors[ion_colors[i]], marker = ".")
            
            ax3.vlines([0,1],[0,0],[1,1],color="w")
            ax3.text(0.05, 0.75, "glysite=%d  %s  %s %d+" %(gpsm.glysite+1, plotmod,
                gpsm.spec, gpsm.charge))
            ax3.get_xaxis().set_ticks([])
            ax3.get_yaxis().set_ticks([])
            
            iwidth = 0
            for igly in range(len(gpsm.glycan)):
                if gpsm.glycan[igly] != 0: iwidth += 1
            LadderStart = 0.025 + (iwidth+1)*2*vfactor
            font_size = 28
            font_name = "Courier New"
            
            if gpsm.peptide != "Z":
                
                ins_rotation = '.'
                ins_no_rot = ' '
                
                colpeptide = " "*len(plotpeptide)
                refpeptide = gpsm.peptide
                glysitepeptide = colpeptide[:gpsm.glysite] + refpeptide[gpsm.glysite] + colpeptide[gpsm.glysite+1:]
                plotpeptide = plotpeptide[:gpsm.glysite] + " " + plotpeptide[gpsm.glysite+1:]
                for modidx, modname in modlist:
                    if modidx == len(gpsm.peptide)+1: modidx -= 2
                    elif modidx > 0: modidx = modidx-1
                    colpeptide = colpeptide[:modidx] + refpeptide[modidx] + colpeptide[modidx+1:]
                    plotpeptide = plotpeptide[:modidx] + " " + plotpeptide[modidx+1:]
                
                ax3.text(LadderStart, 0.25, s=ins_no_rot+plotpeptide+ins_no_rot, fontsize=font_size, fontname=font_name)
                ax3.text(LadderStart, 0.25, s=ins_no_rot+colpeptide+ins_no_rot, fontsize=font_size, fontname=font_name, color="green")
                ax3.text(LadderStart, 0.25, s=ins_no_rot+glysitepeptide+ins_no_rot, fontsize=font_size, fontname=font_name, color="red")
            
            iwidth = 0
            for igly in range(len(gpsm.glycan)):
                if gpsm.glycan[igly] != 0:
                    ax3.plot(0.075 + iwidth*2*vfactor, 0.4, marker=config.used_markers[igly].shape,
                          markerfacecolor=config.used_markers[igly].color, markersize=12, markeredgecolor="black", fillstyle=(config.used_markers[igly].fill), markerfacecoloralt=config.used_markers[igly].alt)
                    ax3.text(0.075 + (iwidth*2+0.5)*vfactor, 0.4, str(gpsm.glycan[igly]),
                            verticalalignment="center",horizontalalignment="left", fontsize = 14)
                    iwidth += 1
            ax1.set_xlim(min_x,max_x)
            ax1.set_ylabel("Intensity")
            ax1.set_xlabel("m/z")
            axLC.set_xlabel("Retention Time (Seconds)")
            axLC.set_ylabel("Intensity")
            
            ax1.xaxis.set_tick_params(color="w")
            
            x1,x2,_y1,_y2 = ax1.axis()
            ax2.axis( (x1,x2,-self.tol, self.tol) )
            
            span = self.tol / 4.
            yidx2 = np.arange(-self.tol, self.tol+span*0.5, span)
            ax2.hlines( yidx2[1:-1], [x1]*(len(yidx2)-2), [x2]*(len(yidx2)-2) ,linestyles = "dashed", colors="gray")
            ax2.hlines( 0, x1, x2 , colors="gray")
            yidx2 = np.arange(-self.tol, self.tol+span, 2*span)
            ax2.set_yticks(yidx2, minor=False)
            ax2.yaxis.set_tick_params(color="w")
            ax2.xaxis.set_tick_params(color="w")
            
            ax2.set_ylabel(r"$\Delta$m (" + self.tol_type + ")")
        return True
    
    def CalcFragmentTol(self, mz):
        if self.tol_type == "Da":
            return self.tol
        else:
            return self.tol * mz / 1e6
    def CalcPPMTol(self, delta, mz):
        if self.tol_type == "Da":
            return delta
        else:
            return delta / mz * 1e6
    #
    def MatchPeak(self, mz_ints, ion_mass):
        isppm = 1 if self.tol_type == "ppm" else 0
        peakidx = spi.PeakIndexing(mz_ints, self.tol, isppm)
        idx, mass_tol, intens = spi.Match(mz_ints, peakidx, ion_mass, self.tol, isppm)
        return idx, mass_tol, intens
    
    def ReadDBSearchRes(self,psm_file):
        self.gpsms.psmlist = {}
        self.gpsms.ReadDBSearchRes(psm_file)
    
    def ReadDenovoRes(self, psm_file):
        self.gpsms.psmlist = {}
        self.gpsms.ReadDenovoRes(psm_file)
        

    def save_plot(self,save_dir):
        start = time.perf_counter() #time.perf_counter
        
        isPlotBack = config.isPlot
        config.isPlot = config.isBatchPlot
        self.output_info = True
        self.outmsg = open(os.path.join(save_dir, self.input_spec+"-glabel.txt"),"w")
        self.outmsg.write("spec\tpeptide\tmodinfo\tglycan(%s)\tformula\tglysite\tcharge\ttheo_ion\tmatched_ion\t%s\n"%(','.join(config.used_glyco),'\t'.join(config.used_glyco)))
        for i, gpsm in enumerate(self.gpsms.psmlist.values()):
            print("[START] %dth GPSM: %s"%(i+1, gpsm.spec))
            _start = time.perf_counter()
            if not gpsm.spec in self.gpsms.psmlist:
                print("no psm of spectrum \"%s\" in result" %spec)
                continue
            
            # for i in range(len(gpsm.peptide)):
                # gpsm.glysite = i
                # plotted = self.SeeOnePlot_new(gpsm)
            plotted = self.SeeOnePlot_new(gpsm)
            if config.isPlot and plotted:
                plt.tight_layout()
                # mng = plt.get_current_fig_manager()
                # mng.window.showMaximized()
                plt.savefig(os.path.join(save_dir, "%s-%s-%s.png"%(gpsm.peptide,config.FormatGlycan(gpsm.glycan), gpsm.spec)),format="png",dpi=config.dpi)
                self.fig.clear()
                plt.close()
            print("[END]   %dth GPSM, %.3f seconds\n"%(i+1, time.perf_counter() - _start))
        self.outmsg.close()
        config.isPlot = isPlotBack
        self.output_info = False
        
        end = time.perf_counter()
        
        print("%d GPSMs, %.3f seconds" %(len(self.gpsms.psmlist), end - start))
        
    def see_oneplot(self, spec):
        if not spec in self.gpsms.psmlist:
            print("no psm of spectrum \"%s\" in result files" %spec)
            return
        gpsm = self.gpsms.psmlist[spec]
        
        # for i in range(len(gpsm.peptide)):
            # gpsm.glysite = i
            # plotted = self.SeeOnePlot_new(gpsm)
        plotted = self.SeeOnePlot_new(gpsm)
        if config.isPlot and plotted:
            plt.tight_layout()
            plt.show()
            self.fig.clear()
            plt.close()

class GUIgLabel(wx.App):
    def OnInit(self):
        self.frame = wx.Frame(parent=None, id=-1,title="gLabel for MS1",
                         pos=(100,100),size=(600,540),
                         style=wx.DEFAULT_FRAME_STYLE,
                         name="frame")
        self.panel = wx.Panel(self.frame,-1)
        
        base_height = 50
        
        self.has_plot = False
        
        self.glabel = Label(tol=10,tol_type="ppm")
        
        label_tol = wx.StaticText(self.panel, -1, "Tolerance:", pos=(40, base_height))
        
        self.tol_type = "ppm"
        self.tolText = wx.TextCtrl(self.panel, -1, pos=(120,base_height),size=(80,-1),style=wx.ALIGN_RIGHT)
        self.tolText.SetValue("10.0")
        self.tolComboBox = wx.ComboBox(self.panel, -1, value="ppm",
                    pos=(220,base_height), choices=["ppm","Da"],
                    style=wx.CB_READONLY)
#         self.tol_update_button = wx.Button(self.panel, -1, 'update', pos=(440,base_height))
#         self.Bind(wx.EVT_BUTTON, self.OnTolUpdate, self.tol_update_button)
        self.Bind(wx.EVT_COMBOBOX, self.OnTolChoose, self.tolComboBox)
        self.Bind(wx.EVT_TEXT, self.OnTolText, self.tolText)
        
        self.plot_glycan = "plot glycan only"
        self.plot_peptide = "plot peptide only"
        self.plot_gp = "plot glycan and peptide"
        
        self.plotComboBox = wx.ComboBox(self.panel, -1, value=self.plot_gp,
                    pos=(360,base_height), choices=[self.plot_gp, self.plot_glycan, self.plot_peptide],
                    style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnPlotChoose, self.plotComboBox)
        
        self.config_button = wx.Button(self.panel, -1, 'Load gLabel config', pos=(440,base_height-40))
        self.Bind(wx.EVT_BUTTON, self.OnConfigButton, self.config_button)
        
        label_mgf = wx.StaticText(self.panel, -1, 'RAW:', pos=(40,base_height+50))
        self.raw_file = wx.TextCtrl(self.panel, -1, pos=(120,base_height+50),size=(300,-1))
        self.raw_file.SetEditable(False)
        self.mgf_button = wx.Button(self.panel, -1, 'browse', pos=(440,base_height+50))
        self.Bind(wx.EVT_BUTTON, self.OnMGFButton, self.mgf_button)
        
        self.type_choose = 0
        
#         self.pGlycoRadio = wx.RadioBox(self.panel, -1, label="pGlyco Type:",
#                     pos=(120,base_height+80), choices=["pGlycoDB","pGlycoDenovo"],
#                     style=wx.RA_SPECIFY_COLS)
#         self.Bind(wx.EVT_RADIOBOX, self.OnTypeChoose, self.pGlycoRadio)

        label_result = wx.StaticText(self.panel, -1, 'pGlycoRes:', pos=(40,base_height+100))
        self.res_file = wx.TextCtrl(self.panel, -1, pos=(120,base_height+100),size=(300,-1))
        self.res_file.SetEditable(False)
        self.res_button = wx.Button(self.panel, -1, 'browse', pos=(440,base_height+100))
        self.Bind(wx.EVT_BUTTON, self.OnResButton, self.res_button)
        
#         self.init_button = wx.Button(self.panel, -1, "Init gLabel",pos=(150,200))
#         self.Bind(wx.EVT_BUTTON, self.OnInitButton, self.init_button)
        label_mz = wx.StaticText(self.panel, -1, 'MaxPlotMZ:', pos=(40,base_height+150))
        self.max_mz = wx.TextCtrl(self.panel, -1, pos=(120,base_height+150),size=(80,-1), style=wx.ALIGN_RIGHT, value="2100.0")
        self.activationComboBox = wx.ComboBox(self.panel, -1, value="HCD",
                    pos=(220,base_height+150), choices=["HCD","ETD","ETHCD"],
                    style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnActivationChoose, self.activationComboBox)
        
        self.ShowMassCheckBox = wx.CheckBox(self.panel, -1, label="show mass",
                    pos=(320,base_height+150+5))
        self.ShowMassCheckBox.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.OnMassCheck, self.ShowMassCheckBox)
        
        label_spec = wx.StaticText(self.panel, -1, 'Spectrum:', pos=(40,base_height+200))
        self.spec_name = wx.TextCtrl(self.panel, -1, pos=(120,base_height+200),size=(300,-1))
        self.spec_button = wx.Button(self.panel, -1, 'show', pos=(440,base_height+200))
        self.Bind(wx.EVT_BUTTON, self.OnSpecButton, self.spec_button)
        
        box = wx.StaticBox(self.panel, -1, 'self defined glycopeptide', pos=(30, base_height+240), size = (510,120))
        label_glycan = wx.StaticText(self.panel, -1, 'Glycan:', pos=(40, base_height+270))
        self.glycan = wx.TextCtrl(self.panel, -1, pos=(120, base_height+270), size=(300,-1))
        label_peptide = wx.StaticText(self.panel, -1, 'Peptide:', pos=(40, base_height+320))
        self.peptide = wx.TextCtrl(self.panel, -1, pos=(120, base_height+320), size=(300,-1))
        self.pep_button = wx.Button(self.panel, -1, 'show this', pos=(440,base_height+320))
        self.Bind(wx.EVT_BUTTON, self.OnPepButton, self.pep_button)
        
        label_batch = wx.StaticText(self.panel, -1, 'BatchOut:', pos=(40,base_height+380))
        self.batch_folder = wx.TextCtrl(self.panel, -1, pos=(120,base_height+380),size=(300,-1))
        self.batch_button = wx.Button(self.panel, -1, 'batch', pos=(440,base_height+380))
        self.Bind(wx.EVT_BUTTON, self.OnBatchBrowseButton, self.batch_button)
        
        self.frame.Show()
        return True
        
    def OnMassCheck(self, event):
        self.glabel.show_mass = self.ShowMassCheckBox.GetValue()
        
    def OnConfigButton(self, event):
        openFile=wx.FileDialog(self.panel, "Open conf file", "", "",
                       "conf files (*.gconf)|*.gconf", 
                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if openFile.ShowModal() == wx.ID_CANCEL:
            return
        config.ReadConfig(openFile.GetPath())
    
    def OnTolChoose(self, event):
        self.glabel.tol_type = self.tolComboBox.GetValue()
    
    def OnPlotChoose(self, event):
        if self.plotComboBox.GetValue() == self.plot_glycan:
            self.glabel.plot_glycan = True
            self.glabel.plot_peptide = False
        elif self.plotComboBox.GetValue() == self.plot_peptide:
            self.glabel.plot_glycan = False
            self.glabel.plot_peptide = True
        else:
            self.glabel.plot_glycan = True
            self.glabel.plot_peptide = True
            
    def OnActivationChoose(self, event):
        config.activation_type = self.activationComboBox.GetValue()
        
    def OnTolText(self, event):
        self.glabel.tol = float(self.tolText.GetValue())
    
    def OnTypeChoose(self, event):
        self.type_choose = self.pGlycoRadio.GetSelection()
    
    def OnBatchBrowseButton(self, event):
        openDir = wx.DirDialog(self.panel, "Choose output folder","",
                                wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if openDir.ShowModal() == wx.ID_CANCEL:
            return
        self.batch_folder.SetValue(openDir.GetPath())
        self.glabel.max_plot_mz = self.GetFloat(self.max_mz.GetValue())
        self.glabel.save_plot(self.batch_folder.GetValue())
        print("****  Finish batch output ****")
    
    def OnMGFButton(self, event):
        openFile=wx.FileDialog(self.panel, "Open MS1 file", "", "",
                       "MS1 files (*.raw)|*.raw", 
                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if openFile.ShowModal() == wx.ID_CANCEL:
            return
        self.raw_file.SetValue(openFile.GetPath())
        self.glabel.ReadRAW(self.raw_file.GetValue())
    
    def GetFloat(self, s, default = 2100.0):
        try:
            return float(s)
        except ValueError:
            return default
            
    def OnPepButton(self, event):
        if not self.has_plot:
            try:
                gpsm = gPSM()
                gpsm.spec = self.spec_name.GetValue()
                if gpsm.spec == "":
                    print("Spectrum is empty")
                    return
                    
                # should be G1 G2 G3..Gn|N1 N2 N3..Nn[|glyfrag=G1 G2 G3..;G1 G2 G3..;]
                # Ti is glyco short name
                glycan = self.glycan.GetValue().split('|')
                if len(glycan) < 2 or glycan[0].strip() == "":
                    print('Glycan should be "G1 G2 G3..Gn|N1 N2 N3..Nn[|glycan fragments in pGlyco format]"')
                    return
                    
                peptide = self.peptide.GetValue().split('|') # should be sequence|glysite[|modification]
                if len(peptide) < 2:
                    print('Peptide should be "peptide|glysite[|modification in pGlyco format]"')
                    return
                
                gpsm.peptide = peptide[0]
                try:
                    gpsm.glysite = int(peptide[1])-1
                except ValueError:
                    print("glysite should be an integer")
                    return
                if len(peptide) > 2:
                    gpsm.modlist = GetModList(peptide[2])
                else:
                    gpsm.modlist = []
                    
                glyco = glycan[0].strip().split(' ')
                glyco = config.UseGlyco(glyco)
                
                gpsm.glycan = tuple([int(g) for g in glycan[1].strip().split(' ')])
                gpsm.glycan_list = []
                if len(glycan) > 2:
                    for gly in glycan[2].strip(';').split(';'):
                        if gly == "": break
                        gpsm.glycan_list.append(tuple([int(g) for g in gly.strip().split(' ')]))
                for g in config.core_frags:
                    if np.all(np.array(g) <= np.array(gpsm.glycan)):
                        gpsm.glycan_list.append(g)
                gpsm.DeleteDuplicated()
            
                self.has_plot = True
                self.glabel.max_plot_mz = self.GetFloat(self.max_mz.GetValue())
                plotted = self.glabel.SeeOnePlot_new(gpsm)
                if config.isPlot and plotted:
                    plt.tight_layout()
                    plt.show()
                config.UseGlyco(glyco)
            except:
                if glyco is not None: config.UseGlyco(glyco)
                self.Except()
            self.has_plot = False
        else:
            print("gLabel is plotting a GPSM, please close that plot!")
    
    def Except(self):
        plt.close()
        print(traceback.format_exc())
        print("An error occurs, please check your input!")
        
    
    def OnSpecButton(self, event):
        if self.spec_name.GetValue() == "": return
        if not self.has_plot:
            self.has_plot = True
            self.glabel.max_plot_mz = self.GetFloat(self.max_mz.GetValue())
            try:
                self.glabel.see_oneplot(self.spec_name.GetValue())
            except:
                self.Except()
            self.has_plot = False
        else:
            print("gLabel is plotting a GPSM, please close that plot!")
    
    def OnResButton(self, event):
        openFile=wx.FileDialog(self.panel, "Open TXT file", "", "",
                       "text files (*.txt)|*.txt",wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if openFile.ShowModal() == wx.ID_CANCEL:
            return
        self.res_file.SetValue(openFile.GetPath())
        if config.ResultHasFragments:
            self.glabel.ReadDenovoRes(self.res_file.GetValue())
        else:
            self.glabel.ReadDBSearchRes(self.res_file.GetValue())
#         if self.type_choose == 1:
#             self.glabel.ReadDenovoRes(self.res_file.GetValue())
#         else:
#             self.glabel.ReadDBSearchRes(self.res_file.GetValue())

if __name__ == "__main__":
    print(" *** gLabel in pGlyco ***")
    app = GUIgLabel()
    app.MainLoop()

        