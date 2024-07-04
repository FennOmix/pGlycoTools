#coding=utf-8
'''
Created on 2013.12.13

@author: dell
'''

import numpy as np
import matplotlib
# matplotlib.use('Qt5Agg')
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import SpecReader
import AAMass
import os
import wx
import copy
import datetime, time, sys
import traceback
import yaml

from alpharaw.match.match_utils import match_closest_peaks
from chemical import replace_element_and_calc_mass

from math import log10

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def load_yaml(yaml):
    with open(yaml) as f:
        return yaml.load(f, Loader=yaml.FullLoader)

aamass = AAMass.AAMass()

def GetRawScanFromSpec(spec):
    items = spec.split(".")
    if len(items) < 6: return ".".join(items[:-1]), items[-1]
    else: return ".".join(items[:-5]), items[-4]

class Marker:
    def __init__(self, glyco = "H", shape = "o", color = "green", marker = "", alt = "white", fill = "full", label = []):
        self.glyco = glyco
        self.shape = shape
        self.color = color
        self.alt = alt
        self.fill = fill
        if not marker: self.oxonium_ions = []
        else:
            self.oxonium_ions = []
            for i in marker.strip().split(";"):
                if i[0].isdigit(): self.oxonium_ions.append(eval(i))
                else:
                    self.oxonium_ions.append(replace_element_and_calc_mass(i, label)[1]+aamass.mass_proton)

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
        self.aaini = ""
        self.dpi = 120
        self.label = []
        self.SetDefault()
        self.settings = {}
        self.ReadConfig()
        
    def SetDefault(self):
        self.used_markers = []
        self.used_markers.append(Marker('H', 'o', 'green', '145.0495347452;163.0600994315;366.1394719645'))
        self.used_markers.append(Marker('N', 's', 'blue', '138.0552587690;168.0655191604;186.0760838467;204.0866485330'))
        self.used_markers.append(Marker('A', 'D', 'purple', '274.0921278414;292.1026925277;657.2348884922'))
        self.used_markers.append(Marker('G', 'D', 'cyan', '290.08704246349997;308.0976071497999;673.2298031143'))
        self.used_markers.append(Marker('F', '^', 'red', '147.0651848094;350.1445573424'))
        self.used_glyco = ['H','N', 'A', 'G', 'F']
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
            kargs['label'] = self.label
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
            elif line.startswith("label"):
                self.label = [item.split("~") for item in line[line.find("=")+1:].strip().split(",")]
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
            elif line.startswith("aaini"):
                self.aaini = line[line.find("=")+1:].strip()
            elif line.startswith("modini"):
                self.modini = line[line.find("=")+1:].strip()
            # elif line.startswith("AA"):
                # items = line[line.find("=")+1:].strip().split(":")
                # aamass.aa_mass_dict[items[0].strip()] = float(items[1].strip())
            elif "shape:" in line and "color:" in line:
                glyco_shape_color(line)
        self.RearrangeGlyco(self.used_glyco)
        if self.glycoini: aamass.__read_glyco_ini__(self.glycoini)
        if self.modini: aamass.__read_mod_ini__(self.modini)
        if self.aaini: aamass.__read_aa_ini__(self.aaini)

config = CGlobalConfig()

fontsize = 12
markersize = 8

hfactor = 0.08
vfactor = 0.025

match_lw = 1
unmatch_lw = 0.5

oxonium_tol = 0.03

#end parameters

#
def CalcPepMass(peptide, modlist):
    if peptide[0].isdigit(): return float(peptide) 
    pep_mass = aamass.mass_H2O
    modmass = [0]*(len(peptide)+2)
    for mod in modlist:
        modmass[mod[0]] = aamass.mod_mass_dict[mod[1]]
    for char in peptide:
        pep_mass += aamass.aa_mass_dict[char]
    for mass in modmass:
        pep_mass += mass
    return pep_mass

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
#
def CalcGlycanMass(glycan_comp):
    glycanmass = 0
    for i in range(len(config.used_glyco)):
        glycanmass += glycan_comp[i]*aamass.glyco_mass_dict[config.used_glyco[i]]
    return glycanmass
#
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
        self.site_group = ""
        self.precursor_mz = 0
        self.ETD_scan = ""
    
    def DeleteDuplicated(self):
        self.glycan_list = list(set(self.glycan_list))

class gPSMList:
    def __init__(self, db_file = config.glycandb):
        self.GlycanCol = 'Glycan'
        self.psmlist = {}
        self.ETDspec_dict = {}
        
    def ReadDenovoRes(self, psm_file):
        self.psmlist = {}
        f = open(psm_file)
        line = f.readline()
        items = [item.strip('"') for item in line[:-1].split("\t")]
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
            
            gpsm.peptide = items[item_idx["Peptide"]]
            gpsm.glysite = int(items[item_idx["GlySite"]]) - 1
            if "Charge" in item_idx: gpsm.charge = int(items[item_idx["Charge"]])
            else: gpsm.charge = int(glyspec.split('.')[-3])
            gpsm.mod = items[item_idx["Mod"]].strip('"')
            gpsm.modlist = GetModList(gpsm.mod)
            
            if "PrecursorMZ" in item_idx: gpsm.precursor_mz = float(items[item_idx["PrecursorMZ"]])
            
            glycos = items[item_idx[self.GlycanCol]].strip().split(" ")
            gpsm.glycan = tuple(int(glyco) for glyco in glycos)
            glycan_list = items[item_idx["GlyFrag"]].strip(';').split(";")
            gpsm.glycan_list = []
            for glycan in glycan_list:
                if len(glycan) < len(config.used_glyco): continue
                glycos = glycan.strip().split(" ")
                iglycan = tuple( int(glycos[i]) for i in range(len(glycos)) )
                gpsm.glycan_list.append(iglycan)
            config.ExtendCoreFrags(gpsm.glycan,gpsm.glycan_list)
            gpsm.DeleteDuplicated()
            
            if "LocalizedSiteGroups" in item_idx:
                gpsm.site_group = items[item_idx["LocalizedSiteGroups"]]
            if 'ETDScan' in item_idx and item_idx["ETDScan"] != '-1':
                gpsm.ETD_scan = items[item_idx["ETDScan"]]
            
            self.psmlist[gpsm.spec] = gpsm
            if 'Scan' in item_idx and 'RawName' in item_idx:
                gpsm = copy.deepcopy(gpsm)
                gpsm.spec = '%s.%s'%(items[item_idx['RawName']],items[item_idx['Scan']])
                self.psmlist[gpsm.spec] = gpsm
            if gpsm.spec != glyspec:
                gpsm = copy.deepcopy(gpsm)
                gpsm.spec = glyspec
                self.psmlist[gpsm.spec] = gpsm
            
            if "ETDScan" in item_idx and item_idx["ETDScan"] != '-1':
                gpsm.ETD_scan = items[item_idx["ETDScan"]]
                
                if 'RawName' in item_idx:
                    gpsm = copy.deepcopy(gpsm)
                    gpsm.spec = '%s.%s'%(items[item_idx['RawName']],items[item_idx['ETDScan']])
                    self.psmlist[gpsm.spec] = gpsm
                
                gpsm = copy.deepcopy(gpsm)
                gpsm.spec = glyspec.replace('.'+items[item_idx['Scan']]+'.'+items[item_idx['Scan']]+'.', '.'+items[item_idx['ETDScan']]+'.'+items[item_idx['ETDScan']]+'.')
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
        self.specfile = None
        self.colory = "orange"
        self.colorb = "cadetblue"
        self.colorz = "magenta"
        self.colorc = "green"
        self.show_mass = False
        self.use_pGlycoSite = False
        
        self.max_plot_mz = 4100.0
        self.save_format = "eps"
        
    def ReadMGF(self, input_file):
        self.input_spec = os.path.split(input_file)[1]
        if self.specfile is not None: self.specfile.close()
        pf2 = input_file[:-3]+'pf2'
        if input_file.endswith('.raw'):
            self.spec_type = 'raw'
            self.specfile = SpecReader.RawFileReader(input_file)
            self.specfile.close = self.specfile.Close
            self.reader = SpecReader.RawReader(input_file)
        elif os.path.isfile(pf2):
            self.spec_type = 'pf2'
            self.specfile = open(pf2, 'rb')
            self.reader = SpecReader.PF2Reader(pf2)
        else:
            self.spec_type = 'mgf'
            self.specfile = open(input_file)
            self.reader = SpecReader.MGFReader(input_file)
        self.reader.ReadAll()
    
    def SeeOnePlot_new(self, gpsm, xmz=None, yint=None):
        print("%s: %s-%s"%(gpsm.spec, gpsm.peptide, config.FormatGlycan(gpsm.glycan)))
        if not self.reader.Has_Spec(gpsm.spec):
            print("no spectrum named \'%s\' in spectrum files" %gpsm.spec)
            return False
        ms2 = self.reader.Get_MS2(gpsm.spec)
        ms2.mz_ints = self.reader.ReadMS2ByIdx(ms2, self.specfile)
        
        if ms2.charge == 0:
            pre_charge = gpsm.charge
        else:
            pre_charge = ms2.charge
        
        ############### init plot ###############
        mz_ints = np.array(ms2.mz_ints)
        if xmz is None:
            xmz = mz_ints[:,0]
            yint = mz_ints[:,1]
        max_inten_real = np.max(yint)
        max_factor = 1.1
        basepeak_after_mz = np.max(yint[xmz >= config.max_oxonium_mz])
        max_inten = basepeak_after_mz * max_factor
        max_plot_mz = min(np.max(xmz),self.max_plot_mz)
        
        yint[yint > max_inten] = max_inten
        
        yint = yint[xmz <= max_plot_mz]
        xmz = xmz[xmz <= max_plot_mz]
        max_plot_mz += 50 # for plot margin
        
        mz_ints = np.append(xmz.reshape(-1,1), yint.reshape(-1,1), axis=1)
        
        margin_Y = 0.6 # 0.9
        margin_by = 0.4 # 0.7
        
        if config.isPlot:
            gs = gridspec.GridSpec(3, 1, height_ratios=[2,8,1])
            self.fig = plt.figure(figsize=(16,10)) #
            ax3 = self.fig.add_subplot(gs[0,0])
            ax1 = self.fig.add_subplot(gs[1,0])
            ax2 = self.fig.add_subplot(gs[2,0])
            ax1.vlines([0,max_plot_mz],[0.2,0],[max_inten*(1+margin_Y*2+margin_by*2),max_inten*(1+margin_Y*2+margin_by*2)],color="w")
            ax1.hlines(0,0,max_plot_mz)
        
            #ax1.set_ylim([0,max_inten])
            ax1.vlines(xmz, [0]*len(yint), yint, color='dimgrey', linewidth=unmatch_lw)
        ############### end init plot ###############
        
        
        ############### glycan and peptide mass ###############
        pepmass = CalcPepMass(gpsm.peptide, gpsm.modlist)
        glycan_mass = CalcGlycanMass(gpsm.glycan)
        gp_mass = pepmass + glycan_mass
        delta = (ms2.pepmass-aamass.mass_proton)*pre_charge-gp_mass
        peplen = len(gpsm.peptide)
        print("precursor mass=%.3f, pepmass=%.3f, glymass=%.3f" %((ms2.pepmass-aamass.mass_proton)*pre_charge, pepmass, glycan_mass))
        ############### end glycan and peptide mass ###############
        
        
        glyfrag_mass = []
        for _gly_frag in gpsm.glycan_list:
            glyfrag_mass.append(CalcGlycanMass(_gly_frag))
        glyfrag_comp = np.array(gpsm.glycan_list, dtype=np.int32)#
        glyfrag_mass = np.array(glyfrag_mass)
        
        self.glycan_0 = np.array([0]*len(config.used_glyco), dtype=np.int32).reshape((1, len(config.used_glyco)))
        
        def add_ions(ion, ion_type, glycan_comp = None):
            if glycan_comp is None: glycan_comp = np.repeat(self.glycan_0, len(ion), axis=0)
            if self.ions is None:
                self.ions = ion
                self.ion_types = ion_type
                self.glycan_comps = glycan_comp
            else:
                self.ions = np.append(self.ions, ion)
                self.ion_types = np.append(self.ion_types, ion_type)
                self.glycan_comps = np.append(self.glycan_comps, glycan_comp, axis = 0)
        
        def add_charged_ions(ion, ion_type, charge, glycan_comp = None):
            if glycan_comp is None: glycan_comp = np.repeat(self.glycan_0, len(ion), axis=0)
            if self.charged_ions is None:
                self.charged_ions = np.array(ion) + aamass.mass_proton
                self.charged_ion_types = ion_type.copy()
                self.charged_glycan_comps = glycan_comp.copy()
                self.charges = np.array([1]*len(ion))
                for i in range(2, charge+1):
                    self.charged_ions = np.append(self.charged_ions, np.array(ion)/i + aamass.mass_proton)
                    self.charged_ion_types = np.append(self.charged_ion_types, ion_type)
                    self.charged_glycan_comps = np.append(self.charged_glycan_comps, glycan_comp, axis = 0)
                    self.charges = np.append(self.charges, np.array([i]*len(ion)))
            else:
                for i in range(1, charge+1):
                    self.charged_ions = np.append(self.charged_ions, np.array(ion)/i + aamass.mass_proton)
                    self.charged_ion_types = np.append(self.charged_ion_types, ion_type)
                    self.charged_glycan_comps = np.append(self.charged_glycan_comps, glycan_comp, axis = 0)
                    self.charges = np.append(self.charges, np.array([i]*len(ion)))
        
        
        by_max_charge = 1 if pre_charge <= 2 else 2
        cz_max_charge = 1 if pre_charge <= 2 else 2
        Y_max_charge = pre_charge - 1
        M_max_charge = pre_charge
        
        self.ions = None
        self.charged_ions = None
        
        ############### Y ions ###############
        if len(glyfrag_mass) > 0:
            add_ions(glyfrag_mass + pepmass, np.array(['Y']*len(glyfrag_mass)), glyfrag_comp)
            add_charged_ions(self.ions, self.ion_types, Y_max_charge, self.glycan_comps)
        else:
            #empty Y ions
            add_ions(np.array([0]), np.array(['']), self.glycan_0)
        if gpsm.peptide != "Z":
            Y0_ion = [pepmass]
            Y0_ion_type = ['Y0']
            add_ions(Y0_ion, Y0_ion_type)
            add_charged_ions(Y0_ion, Y0_ion_type, Y_max_charge)
        
        M_ion = [gp_mass]
        M_ion_type = ['M']
        add_ions(M_ion, M_ion_type)
        add_charged_ions(M_ion, M_ion_type, M_max_charge)
        ############### end Y ions ###############
        
        if gpsm.peptide == "Z":
            ############### B ions ###############
            B_ions = glyfrag_mass
            B_ion_types = np.array(['B']*len(glyfrag_mass))
            add_ions(B_ions, B_ion_types, glyfrag_comp)
            add_charged_ions(B_ions, B_ion_types, Y_max_charge, glyfrag_comp)
            
            M_ion = [glycan_mass]
            M_ion_type = ['M-H2O']
            add_ions(M_ion, M_ion_type)
            add_charged_ions(M_ion, M_ion_type, M_max_charge)
        ############### end Y ions ###############
        
        ############### base b/y, c/z mass ###############
        if config.activation_type.upper() == "ETHCD" or config.activation_type.upper() == "HCDPDETXXD":
            has_by = True
            has_cz = True
        elif config.activation_type.upper() == "ETD":
            has_by = False
            has_cz = True
        else:
            has_by = True
            has_cz = False
        if gpsm.peptide == "Z":
            has_by = False
            has_cz = False
        
        modmass = [0]*(len(gpsm.peptide)+2)
        for mod in gpsm.modlist:
            modmass[mod[0]] += aamass.mod_mass_dict[mod[1]]
        
        def get_b_y_ions():
            bion = []
            yion = []
            mass_nterm = modmass[0]
            i = 1
            for char in gpsm.peptide[:-1]:
                mass_nterm += aamass.aa_mass_dict[char] + modmass[i]
                i += 1
                bion.append(mass_nterm)
                yion.append(pepmass - mass_nterm)
            return np.array(bion), np.array(yion)
        bion, yion = get_b_y_ions()
            
        max_glyco_in_by = 1
        if has_by:
            # self.glycan_comps = np.repeat(self.glycan_0, len(bion), axis=0)
            b_ion_type = ["b%d"%(i+1) for i in range(len(bion))]
            add_ions(bion, b_ion_type)
            add_charged_ions(bion, b_ion_type, by_max_charge)
            y_ion_type = ["y%d"%(peplen-i-1) for i in range(len(yion))]
            add_ions(yion, y_ion_type)
            add_charged_ions(yion, y_ion_type, by_max_charge)
        
        
            if config.glycan_type == 'N':
                ############### b/y + $ (cxr = cross-ring) ###############
                _cxr_mass = 83.03711
                bion_cxr = bion[gpsm.glysite:] + _cxr_mass
                bion_cxr_type = ["b$%d"%(i+1) for i in range(gpsm.glysite,len(bion))]
                add_ions(bion_cxr, bion_cxr_type)
                add_charged_ions(bion_cxr, bion_cxr_type, by_max_charge)
                
                yion_cxr = yion[:gpsm.glysite] + _cxr_mass
                yion_cxr_type = ["y$%d"%(peplen-i-1) for i in range(gpsm.glysite)]
                add_ions(yion_cxr, yion_cxr_type)
                add_charged_ions(yion_cxr, yion_cxr_type, by_max_charge)
                
                add_ions([pepmass + _cxr_mass], ['Y$'])
                add_charged_ions([pepmass + _cxr_mass], ['Y$'], Y_max_charge)
                ############### end b/y + $ (cxr = cross-ring) ###############
        
            ############### b/y + glyco unit ###############
            for _comp,_mass in zip(glyfrag_comp, glyfrag_mass):
                if np.sum(_comp) <= max_glyco_in_by:
                    bion_add_glyco = bion[gpsm.glysite:] + _mass
                    bion_add_glyco_type = ["b%d"%(i+1) for i in range(gpsm.glysite,len(bion))]
                    bion_add_glyco_comp = np.repeat(_comp.reshape(1,-1), len(bion_add_glyco), axis=0)
                    add_ions(bion_add_glyco, bion_add_glyco_type, bion_add_glyco_comp)
                    add_charged_ions(bion_add_glyco, bion_add_glyco_type, by_max_charge, bion_add_glyco_comp)
                    
                    yion_add_glyco = yion[:gpsm.glysite] + _mass
                    yion_add_glyco_type = ["y%d"%(peplen-i-1) for i in range(gpsm.glysite)]
                    yion_add_glyco_comp = np.repeat(_comp.reshape(1,-1), len(yion_add_glyco), axis=0)
                    add_ions(yion_add_glyco, yion_add_glyco_type, yion_add_glyco_comp)
                    add_charged_ions(yion_add_glyco, yion_add_glyco_type, by_max_charge, yion_add_glyco_comp)
            ############### end b/y + glyco unit ###############
            
        if self.use_pGlycoSite and gpsm.site_group:
            if "{" in gpsm.site_group: site_groups = gpsm.site_group.strip("}").split("}")
            else: site_groups = gpsm.site_group.strip(";").split(";")
            sites = []
            site_glycans = []
            site_probs = []
            site_glymasses = [0]*len(gpsm.peptide)
            for site_group in site_groups:
                if "{" in gpsm.site_group: items = site_group.strip("{").split(",")
                else: items = site_group.split(",")
                sites.append((int(items[0][1:])-1, int(items[1][1:])-1))
                site_glycans.append([int(glyco) for glyco in items[2][1:-1].split(' ')])
                site_probs.append(items[3])
                site_glymasses[sites[-1][0]] = CalcGlycanMass(site_glycans[-1])
            site_glymasses = np.cumsum(site_glymasses)
            
        if has_cz and self.use_pGlycoSite and gpsm.site_group:
            cion = bion + aamass.mass_NH3 + site_glymasses[:-1]
            cion_type = ["c%d"%(i+1) for i in range(len(cion))]
            add_ions(cion, cion_type)
            add_charged_ions(cion, cion_type, cz_max_charge)
            
            zion = pepmass + glycan_mass + aamass.mass_H - cion
            zion_type = ["z%d"%(peplen-i-1) for i in range(len(zion))]
            add_ions(zion, zion_type)
            add_charged_ions(zion, zion_type, cz_max_charge)
            
            # c_minus_H_ion = cion - aamass.mass_H
            # cH_ion_type = ["c%d-H"%(i+1) for i in range(len(cion))]
            # add_ions(c_minus_H_ion, cH_ion_type)
            # add_charged_ions(c_minus_H_ion, cH_ion_type, cz_max_charge)
            
            z_plus_H_ion = zion + aamass.mass_H
            zH_ion_type = ["z%d+H"%(peplen-i-1) for i in range(len(zion))]
            add_ions(z_plus_H_ion, zH_ion_type)
            add_charged_ions(z_plus_H_ion, zH_ion_type, cz_max_charge)
        elif has_cz:
            cion = bion + aamass.mass_NH3
            cion[gpsm.glysite:] += glycan_mass
            cion_type = ["c%d"%(i+1) for i in range(len(cion))]
            add_ions(cion, cion_type)
            add_charged_ions(cion, cion_type, cz_max_charge)
            
            zion = yion - aamass.mass_NH3 + aamass.mass_H
            zion[:gpsm.glysite] += glycan_mass
            zion_type = ["z%d"%(peplen-i-1) for i in range(len(zion))]
            add_ions(zion, zion_type)
            add_charged_ions(zion, zion_type, cz_max_charge)
            
            # c_minus_H_ion = cion - aamass.mass_H
            # cH_ion_type = ["c%d-H"%(i+1) for i in range(len(cion))]
            # add_ions(c_minus_H_ion, cH_ion_type)
            # add_charged_ions(c_minus_H_ion, cH_ion_type, cz_max_charge)
            
            z_plus_H_ion = zion + aamass.mass_H
            zH_ion_type = ["z%d+H"%(peplen-i-1) for i in range(len(zion))]
            add_ions(z_plus_H_ion, zH_ion_type)
            add_charged_ions(z_plus_H_ion, zH_ion_type, cz_max_charge)
        ############### base b/y, c/z mass ###############
        
        
        if self.output_info:
            #spec\tpeptide\tmodinfo\tglycan(%s)\tformula\tglysite\tcharge\ttheo_ion\tmatched_ion
            self.outmsg.write("%s\t%s\t%s" %(gpsm.spec, gpsm.peptide, gpsm.mod))
            self.outmsg.write("\t%s"%",".join([str(g) for g in gpsm.glycan]))
            self.outmsg.write("\t%s"%config.FormatGlycan(gpsm.glycan))
            if self.use_pGlycoSite and gpsm.site_group:
                self.outmsg.write("\t%s"%gpsm.site_group)
            else:
                self.outmsg.write("\t%d"%(gpsm.glysite+1))
            self.outmsg.write("\t%d"%(gpsm.charge))
            output_strs = []
            for i in range(len(self.ions)):
                if self.ion_types[i] != "Y" and (sum(self.glycan_comps[i]) > max_glyco_in_by): continue
                gly_format = config.FormatGlycan(self.glycan_comps[i])
                if gly_format != "": gly_format = '-' + gly_format
                output_strs.append("%s%s"%(self.ion_types[i],gly_format))
            self.outmsg.write("\t%s;"%(";".join(output_strs)))
        
        clip = self.charged_ions <= (self.max_plot_mz + 1)
        
        ions = self.charged_ions[clip]
        glycan_comps = self.charged_glycan_comps[clip] 
        ion_types = self.charged_ion_types[clip]
        charges = self.charges[clip]
        
        idx, mass_tol = self.MatchPeak(mz_ints, ions)
        idx = np.array(idx)
        mass_tol = np.array(mass_tol)
        
        ions = ions[idx != -1]
        glycan_comps = glycan_comps[idx != -1]
        ion_types = ion_types[idx != -1]
        charges = charges[idx != -1]
        mass_tol = mass_tol[idx != -1]
        idx = idx[idx != -1]
        
        # output the matched ion type
        if self.output_info:
            output_strs = []
            for i in range(len(ions)):
                if ion_types[i] != "Y" and (sum(glycan_comps[i]) > max_glyco_in_by): continue
                gly_format = config.FormatGlycan(glycan_comps[i])
                if gly_format != "": gly_format = '-' + gly_format
                output_strs.append("%s%s+%d=%.3f,%.1f" %(ion_types[i],gly_format, charges[i], xmz[idx[i]], yint[idx[i]]))
            self.outmsg.write("\t%s"%(";".join(output_strs)))
        
        def _match_oxonium(num, markers, color):
            tol = oxonium_tol
            output_strs = []
            matched_oxonium = []
            if num > 0:
                for mass in markers:
                    matched_idx = -1
                    for i in range(len(mz_ints)):
                        if mz_ints[i,0] > mass+tol: break
                        elif abs(mass - mz_ints[i,0]) <= tol:
                            if matched_idx == -1 or mz_ints[matched_idx,1] < mz_ints[i,1]:
                                matched_idx = i
                    i = matched_idx
                    if i != -1:
                        output_strs.append("%.2f=%.1f"%(mz_ints[i,0], mz_ints[i,1]))
                        matched_oxonium.append((mz_ints[i,0], mz_ints[i,1], color))
            if self.output_info:
                self.outmsg.write("\t%s"%(";".join(output_strs)))
            return matched_oxonium
            
        matched_oxoniums = []
        for i in range(len(config.used_markers)):
            matched_oxoniums.append(_match_oxonium(gpsm.glycan[i], config.used_markers[i].oxonium_ions, config.used_markers[i].color))
        
        if self.output_info:
            self.outmsg.write("\n")
            self.outmsg.flush()

        if config.isPlot:
            ########## plot
            order = np.argsort(idx)
            mass_tol = mass_tol[order]
            ions = ions[order]
            glycan_comps = glycan_comps[order]
            ion_types = ion_types[order]
            charges = charges[order]
            idx = idx[order]
            peptide = gpsm.peptide
            plotpeptide = peptide
            peplen = len(peptide)
            plotmod = "Mod: "
            for modidx, modname in gpsm.modlist:
                if modidx == 0: plotmod += "NTerm" + aamass.mod_to_mass[modname] + ";"
                elif modidx == peplen+1: plotmod += "CTerm" + aamass.mod_to_mass[modname] + ";"
                else: plotmod += gpsm.peptide[modidx-1] + str(modidx) + aamass.mod_to_mass[modname] + ";"
                modidx = modidx - peplen
                if modidx >= 0:
                    peptide = peptide + aamass.mod_to_mass[modname]
                else:
                    peptide = peptide[-len(peptide):modidx] + aamass.mod_to_mass[modname] + peptide[modidx:]
            if len(gpsm.modlist) == 0:
                plotmod = "noPepMod"
            
            xmztol = xmz
            xmz = xmz[idx]
            yint = yint[idx]
            
            #
            b_ion_set = set()
            y_ion_set = set()
            if self.plot_glycan and self.plot_peptide:
                baseheights = [max_inten, max_inten*(1+margin_Y)]
                bybaseheights = [max_inten*(1+margin_Y*2+margin_by), max_inten*(1+margin_Y*2)]
            elif self.plot_glycan:
                baseheights = [max_inten*(1+margin_Y*0.1), max_inten*(1+margin_Y*1.3)]
            else:
                baseheights = [max_inten, max_inten*(1+margin_Y*0.5)]
                bybaseheights = baseheights
            
            count = 0
            countby = 0
            text_offset = 0
            text_to_mark = .95
            
            # plot oxoniums
            for oxonium in matched_oxoniums:
                for mz,inten,color in oxonium:
                    ax1.text(mz*(1+text_offset),inten+hfactor*max_inten,
                        "%.2f" %mz, rotation = 90, color = color,
                        horizontalalignment="center",verticalalignment="bottom")
                    ax1.plot([mz, mz], [0, inten], color = color, linewidth=match_lw)
            
            cz_set = set()
            for i in range(len(ions)):
                height = 0.1
                if ion_types[i][0] == "Y" or ion_types[i][0] == "B":
                    if not self.plot_glycan: continue
                    elif ion_types[i][0] == "B" and xmz[i] < config.max_oxonium_mz: continue
                    baseheight = baseheights[count%2]
                    count += 1
                    ax1.plot([xmz[i],xmz[i]],[yint[i], baseheight], dashes=[4, 4], color="gray", linewidth=0.2)
                    if not config.glyco_as_text:
                        for igly in range(len(config.used_markers)):
                            if glycan_comps[i][igly] > 0:
                                ax1.plot(xmz[i],baseheight+height*hfactor*max_inten,
                                        marker=config.used_markers[igly].shape, markerfacecolor=config.used_markers[igly].color, markersize=markersize, markeredgecolor="black", fillstyle=fill_rotation(config.used_markers[igly].fill), markerfacecoloralt=config.used_markers[igly].alt)
                                ax1.text(xmz[i],baseheight+((height+text_to_mark)*hfactor)*max_inten,
                                         str(glycan_comps[i][igly]),rotation=90, fontsize=fontsize,
                                         horizontalalignment="center",verticalalignment="bottom")
            
                                height += int(log10(glycan_comps[i][igly])+1)*1.3 + 1.5
                        ax1.text(xmz[i],baseheight+((height)*hfactor)*max_inten,
                                 "%s /%d+%s"%(ion_types[i],charges[i]," %.2f"%xmz[i] if self.show_mass else ""), rotation=90, fontsize=fontsize,
                                 horizontalalignment="center",verticalalignment="bottom")
                                 
                    else:
                        gly_str = config.short_format(glycan_comps[i])
                        if gly_str != "": gly_str = "--" + gly_str
                        gly_str = "%s%s /%d+" %(ion_types[i], gly_str, charges[i])
                        ax1.text(xmz[i], baseheight+height*hfactor*max_inten, gly_str, rotation=90, fontsize=fontsize, horizontalalignment="center",verticalalignment="bottom")
                    
                    ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color="red", linewidth=match_lw)
                    ax2.plot(xmz[i], mass_tol[i], color="red", marker = ".")
                elif ion_types[i][0] == "M":
                    baseheight = baseheights[count % 2]
                    count += 1
                    ax1.plot([xmz[i],xmz[i]],[yint[i], baseheight], dashes=[4, 4], color="gray", linewidth=0.2)
                    ax1.text(xmz[i]*(1+text_offset),baseheight+((height)*hfactor)*max_inten,
                             "%s /%d+%s"%(ion_types[i],charges[i]," %.2f"%xmz[i] if self.show_mass else ""),rotation=90, fontsize=fontsize,
                             horizontalalignment="center",verticalalignment="bottom")
                    ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color="blue", linewidth=match_lw)
                    ax2.plot(xmz[i], mass_tol[i], color="blue", marker = ".")
                else: #b, y ions
                    #Just plot Y1 fragment b,y ions
                    if not self.plot_peptide: continue
                    if sum(glycan_comps[i]) > max_glyco_in_by: continue
                    baseheight = bybaseheights[count % 2]
                    
                    if ion_types[i].startswith('z'):
                        if ion_types[i].endswith('H'):
                            ch_ion = ion_types[i][1:-2] + str(charges[i])
                        else:
                            ch_ion = ion_types[i][1:] + str(charges[i])
                        if ch_ion in cz_set: 
                            continue
                        else:
                            cz_set.add(ch_ion)
                    count += 1
                    ax1.plot([xmz[i],xmz[i]],[yint[i], baseheight], dashes=[4, 4], color="gray", linewidth=0.2)

                    if not config.glyco_as_text:
                        for igly in range(len(config.used_markers)):
                            if glycan_comps[i][igly] > 0:
                                ax1.plot(xmz[i],baseheight+height*hfactor*max_inten,
                                            marker=config.used_markers[igly].shape, markerfacecolor=config.used_markers[igly].color, markersize=markersize, markeredgecolor="black", fillstyle=fill_rotation(config.used_markers[igly].fill), markerfacecoloralt=config.used_markers[igly].alt)
                                ax1.text(xmz[i]*(1+text_offset),baseheight+((height+text_to_mark)*hfactor)*max_inten,
                                         str(glycan_comps[i][igly]),rotation=90, fontsize=fontsize,
                                         horizontalalignment="center",verticalalignment="bottom")
                                         #horizontalalignment="center",verticalalignment="top")
            
                                height += 2.5
                        ax1.text(xmz[i]*(1+text_offset),baseheight+((height)*hfactor)*max_inten, "%s /%d+%s"%(ion_types[i],charges[i]," %.2f"%xmz[i] if self.show_mass else ""),rotation=90, fontsize=fontsize, horizontalalignment="center",verticalalignment="bottom")
                    else:
                        gly_str = config.short_format(glycan_comps[i])
                        if gly_str != "": gly_str = "--" + gly_str
                        gly_str = "%s%s /%d+" %(ion_types[i], gly_str, charges[i])
                        ax1.text(xmz[i], baseheight+height*hfactor*max_inten, gly_str, rotation=90, fontsize=fontsize, horizontalalignment="center",verticalalignment="bottom")
                    
                    if ion_types[i][:2] == "b$":
                        b_ion_set.add(int(ion_types[i][2:]))
                        ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=self.colorb, linewidth=match_lw)
                        ax2.plot(xmz[i], mass_tol[i], color=self.colorb, marker = ".")
                    elif ion_types[i][0] == "b":
                        b_ion_set.add(int(ion_types[i][1:]))
                        ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=self.colorb, linewidth=match_lw)
                        ax2.plot(xmz[i], mass_tol[i], color=self.colorb, marker = ".")
                    elif ion_types[i][:2] == "y$":
                        y_ion_set.add(int(ion_types[i][2:]))
                        ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=self.colory, linewidth=match_lw)
                        ax2.plot(xmz[i], mass_tol[i], color=self.colory, marker = ".")
                    elif ion_types[i][0] == "y":
                        y_ion_set.add(int(ion_types[i][1:]))
                        ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=self.colory, linewidth=match_lw)
                        ax2.plot(xmz[i], mass_tol[i], color=self.colory, marker = ".")
                    elif ion_types[i][0] == "c":
                        if ion_types[i][-1] == "H":
                            b_ion_set.add(int(ion_types[i][1:-2]))
                        else:
                            b_ion_set.add(int(ion_types[i][1:]))
                        ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=self.colorc, linewidth=match_lw)
                        ax2.plot(xmz[i], mass_tol[i], color=self.colorc, marker = ".")
                    elif ion_types[i][0] == "z":
                        if ion_types[i][-1] == "H":
                            y_ion_set.add(int(ion_types[i][1:-2]))
                        else:
                            y_ion_set.add(int(ion_types[i][1:]))
                        ax1.plot([xmz[i],xmz[i]],[0, yint[i]], color=self.colorz, linewidth=match_lw)
                        ax2.plot(xmz[i], mass_tol[i], color=self.colorz, marker = ".")
            
            ax3.vlines([0,1],[0,0],[1,1],color="w")
            delta_ppm = delta/gp_mass*1000000
            def site_group_str(site_glycans, sites):
                _ret = []
                for (i,j),gly,prob in zip(sites, site_glycans, site_probs):
                    if i != j:
                        _ret.append('{}{}-{}{}:{}@{}'.format(gpsm.peptide[i], i+1, gpsm.peptide[j], j+1, config.short_format(gly), prob[2:4]))
                    else:
                        _ret.append('{}{}:{}@{}'.format(gpsm.peptide[i],i+1,config.short_format(gly), prob[2:4]))
                return ",".join(_ret)
            if self.use_pGlycoSite and gpsm.site_group:
                glysite = site_group_str(site_glycans, sites)
            else:
                glysite = str(gpsm.glysite+1)
                sites = [(gpsm.glysite, gpsm.glysite)]
            ax3.text(-0.03, 0.75, "Site=%s %s\n%s %d+" %(glysite, plotmod,
                gpsm.spec, pre_charge) + r" $\Delta$m=%.2f ppm, %.2f Th" %(delta_ppm, delta/pre_charge))
            ax3.get_xaxis().set_ticks([])
            ax3.get_yaxis().set_ticks([])
            
            iwidth = 0
            for igly in range(len(gpsm.glycan)):
                if gpsm.glycan[igly] != 0: iwidth += 1
            LadderStart = (iwidth+1)*2*vfactor
            font_size = 28
            font_name = "Courier New"
            
            def get_unicode(code):
                return code.decode('unicode-escape')
                
            y_ladder = ""
            b_ladder = ""
            def ladder_idx(i):
                return r"$_{_{_{%d}}}$"%i
            
                    
            def plot_ladder_idx(ax, x, y, idx, fontsize, fontname, color):
                ax.text(x, y, idx, fontsize=fontsize,fontname=fontname,color=color)
            
            if gpsm.peptide != "Z":
                subscript=3
                y_idx = ""
                b_idx = ""
                for i in range(1, len(gpsm.peptide)):
                    if i in y_ion_set:
                        y_ladder = get_unicode(b"\u2310") + y_ladder
                        y_idx = "%d"%i+" "*(subscript-1-int(log10(i))) + y_idx
                    else:
                        y_ladder = " " + y_ladder
                        y_idx = " "*subscript + y_idx
                    if i in b_ion_set:
                        b_ladder += get_unicode(b"\u2310")
                        b_idx += " "*(subscript-1-int(log10(i))) + "%d"%i
                    else:
                        b_ladder += " "
                        b_idx += " "*subscript
                y_ladder = " " + y_ladder
                b_ladder = b_ladder + " "
                b_ladder = b_ladder[::-1]#for reverse plot
                
                ins_rotation = '.'
                ins_no_rot = ' '
                
                colpeptide = " "*len(plotpeptide)
                refpeptide = gpsm.peptide
                glysitepeptide = colpeptide
                sitegrouppeptide = colpeptide
                for site,site2 in sites:
                    if site == site2:
                        glysitepeptide = glysitepeptide[:site] + refpeptide[site] + glysitepeptide[site+1:]
                    else:
                        sitegrouppeptide = sitegrouppeptide[:site] + refpeptide[site:site2+1] + sitegrouppeptide[site2+1:]
                    plotpeptide = plotpeptide[:site] + " "*(site2-site+1) + plotpeptide[site2+1:]
                for modidx, modname in gpsm.modlist:
                    if modidx == len(gpsm.peptide)+1: modidx -= 2
                    elif modidx > 0: modidx = modidx-1
                    colpeptide = colpeptide[:modidx] + refpeptide[modidx] + colpeptide[modidx+1:]
                    plotpeptide = plotpeptide[:modidx] + " " + plotpeptide[modidx+1:]
                    
                if has_by:
                    color_cterm = self.colory
                    color_nterm = self.colorb
                else:
                    color_cterm = self.colorz
                    color_nterm = self.colorc
                
                ax3.text(LadderStart, 0.25, s=ins_no_rot+plotpeptide+ins_no_rot, fontsize=font_size, fontname=font_name)
                ax3.text(LadderStart, 0.25, s=ins_no_rot+colpeptide+ins_no_rot, fontsize=font_size, fontname=font_name, color="green")
                ax3.text(LadderStart, 0.25, s=ins_no_rot+glysitepeptide+ins_no_rot, fontsize=font_size, fontname=font_name, color="red")
                ax3.text(LadderStart, 0.25, s=ins_no_rot+sitegrouppeptide+ins_no_rot, fontsize=font_size, fontname=font_name, color="darksalmon")
                ax3.text(LadderStart, 0.48, s=ins_no_rot+y_ladder+ins_no_rot, fontsize=font_size, fontname=font_name,color=color_cterm)
                # t = ax3.text(LadderStart, 0.48, s=ins_no_rot+y_ladder+ins_no_rot, fontsize=font_size, fontname=font_name,color=self.colory)
                ax3.text(LadderStart, 0.06, s=ins_rotation+b_ladder+ins_rotation, fontsize=font_size, fontname=font_name,color=color_nterm, rotation=180)
                
                plot_ladder_idx(ax3, LadderStart, 0.60, " "*(2*subscript+1) + y_idx, font_size/subscript, font_name, color_cterm)
                plot_ladder_idx(ax3, LadderStart+0.003, 0.01, " "*(1*subscript-1) + b_idx, font_size/subscript, font_name, color_nterm)
                
                ax3.text(LadderStart, 0.1, s=get_unicode(b'\u25A0')+' '*len(plotpeptide)+get_unicode(b'\u25A0'), fontsize=font_size, fontname=font_name,color="white")
            
            iwidth = 0
            for igly in range(len(gpsm.glycan)):
                if gpsm.glycan[igly] != 0:
                    ax3.plot(0 + iwidth*2*vfactor, 0.4, marker=config.used_markers[igly].shape,
                          markerfacecolor=config.used_markers[igly].color, markersize=12, markeredgecolor="black", fillstyle=(config.used_markers[igly].fill), markerfacecoloralt=config.used_markers[igly].alt)
                    ax3.text(0 + (iwidth*2+0.5)*vfactor, 0.4, str(gpsm.glycan[igly]),
                            verticalalignment="center",horizontalalignment="left", fontsize = 14)
                    iwidth += 1
            
            ax1.set_xlim(0,max_plot_mz)
            if self.plot_glycan and self.plot_peptide:
                ymax = max_inten*(1 + margin_Y*2.2 + margin_by*2.2)
                ax1.set_ylim(bottom=0,top=ymax)
            elif self.plot_glycan:
                ymax = max_inten*(1 + margin_Y*2.5)
                ax1.set_ylim(bottom=0,top=ymax)
            else:
                ymax = max_inten*(1 + margin_Y)
                ax1.set_ylim(bottom=0,top=ymax)
                
            
            ax1.xaxis.set_tick_params(color="w")
            ax1.text(10,ymax+max_inten*0.05, 'x%.1e'%max_inten_real)
            
            if ymax > max_inten_real:
                step = max_inten_real / 5.0
                ticks = np.arange(step, ymax, step)
                labels = np.around(ticks / max_inten_real * 100, -1)
            else:
                step = ymax / 6.0
                ticks = np.arange(step, ymax, step)
                labels = ticks / max_inten_real * 100
            labels = ['%d%%'%int(i) for i in labels]
            ax1.set_yticks(ticks)
            ax1.set_yticklabels(labels)
            
            
            x1,x2,_y1,_y2 = ax1.axis()
            ax2.axis( (x1,x2,-self.tol, self.tol) )
            
            span = self.tol / 4.
            yidx2 = np.arange(-self.tol, self.tol+span*0.5, span)
            ax2.hlines( yidx2[1:-1], [x1]*(len(yidx2)-2), [x2]*(len(yidx2)-2) ,linestyles = "dashed", colors="gray", linewidth=0.2)
            ax2.hlines( 0, x1, x2 , colors="gray")
            yidx2 = np.arange(-self.tol, self.tol+span, 2*span)
            ax2.set_yticks(yidx2, minor=False)
            ax2.yaxis.set_tick_params(color="w")
            ax2.xaxis.set_tick_params(color="w")
            
            ax2.set_xlabel("m/z")
            ax2.set_ylabel(r"$\Delta$m (" + self.tol_type + ")")
            ax1.set_ylabel("Relative Intensity")
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
        if self.tol_type == "ppm":
            ion_tols = ion_mass*self.tol*1e-6
        else:
            ion_tols = np.ones_like(ion_tols)*self.tol
        idx = match_closest_peaks(mz_ints[:,0], mz_ints[:,1], ion_mass, ion_tols)
        matched_masses = mz_ints[idx,0]
        mass_tol = matched_masses-ion_mass
        if self.tol_type == "ppm":
            mass_tol = mass_tol*1e6/matched_masses
        mass_tol[idx==-1] = 2*self.tol
        return idx, mass_tol
    
    def ReadDBSearchRes(self,psm_file):
        self.gpsms.psmlist = {}
        self.gpsms.ReadDBSearchRes(psm_file)
    
    def ReadDenovoRes(self, psm_file):
        self.gpsms.psmlist = {}
        self.gpsms.ReadDenovoRes(psm_file)
        

    def save_plot(self,save_dir):
        start = time.perf_counter() #time.perf_counter
        
        if not os.path.isdir(save_dir):
            try:
                os.makedirs(save_dir)
            except:
                pass
        
        isPlotBack = config.isPlot
        config.isPlot = config.isBatchPlot
        self.output_info = True
        self.outmsg = open(os.path.join(save_dir, self.input_spec+"-glabel.txt"),"w")
        self.outmsg.write("spec\tpeptide\tmodinfo\tglycan(%s)\tformula\tglysite\tcharge\ttheo_ion\tmatched_ion\t%s\n"%(','.join(config.used_glyco),'\t'.join(config.used_glyco)))
        plotted_set = set()
        for i, gpsm in enumerate(self.gpsms.psmlist.values()):
            if not gpsm.spec.endswith('.dta'): continue
            print("[START] %dth GPSM: %s"%(i+1, gpsm.spec))
            _start = time.perf_counter()
            if not gpsm.spec in self.gpsms.psmlist:
                print("no psm of spectrum \"%s\" in result" %spec)
                continue
            if not self.reader.Has_Spec(gpsm.spec):
                gpsm = self.GetETD_GPSM(gpsm.spec, gpsm)
            
            if gpsm.spec in plotted_set: continue
            plotted_set.add(gpsm.spec)
            # for i in range(len(gpsm.peptide)):
                # gpsm.glysite = i
                # plotted = self.SeeOnePlot_new(gpsm)
            plotted = self.SeeOnePlot_new(gpsm)
            if config.isPlot and plotted:
                plt.tight_layout()
                # mng = plt.get_current_fig_manager()
                # mng.window.showMaximized()
                plt.savefig(os.path.join(save_dir, "%s-%s-%s.%s"%(gpsm.peptide,config.FormatGlycan(gpsm.glycan), gpsm.spec,self.save_format)),format=self.save_format,dpi=config.dpi)
                self.fig.clear()
                plt.close()
            print("[END]   %dth GPSM, %.3f seconds\n"%(i+1, time.perf_counter() - _start))
        self.outmsg.close()
        config.isPlot = isPlotBack
        self.output_info = False
        
        end = time.perf_counter()
        
        print("%d GPSMs, %.3f seconds" %(len(self.gpsms.psmlist), end - start))
        
    def GetETD_GPSM(self, spec, gpsm):
        raw, scan = GetRawScanFromSpec(spec)
        if gpsm.ETD_scan and gpsm.ETD_scan != scan and gpsm.ETD_scan != "-1":
            if ".{scan}.{scan}.".format(scan=scan) in spec:
                spec = spec.replace(".{scan}.{scan}.".format(scan=scan), ".{scan}.{scan}.".format(scan=gpsm.ETD_scan))
            else:
                spec = spec.replace("."+scan, "."+gpsm.ETD_scan)
            if spec in self.gpsms.psmlist:
                print("Change scan={} to ETD scan={} for pGlycoSite".format(scan, gpsm.ETD_scan))
                return self.gpsms.psmlist[spec]
            else:
                print("Change scan={} to ETD scan={}, ETD scan is not in results".format(scan, gpsm.ETD_scan))
                return gpsm
        else:
            return gpsm
        
    def see_oneplot(self, spec):
        if not spec in self.gpsms.psmlist:
            print("no psm of spectrum \"%s\" in result files" %spec)
            return
        gpsm = self.gpsms.psmlist[spec]
        if not self.reader.Has_Spec(spec):
            gpsm = self.GetETD_GPSM(spec, gpsm)
        # for i in range(len(gpsm.peptide)):
            # gpsm.glysite = i
            # plotted = self.SeeOnePlot_new(gpsm)
        plotted = self.SeeOnePlot_new(gpsm)
        if config.isPlot and plotted:
            plt.tight_layout()
            plt.show()
            self.fig.clear()
            plt.close()
            
class GUIgLabel(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((600, 540))
        self.SetTitle("gLabel for glycopeptide")
        
        self.panel = wx.Panel(self,-1)
        
        base_height = 50
        
        self.has_plot = False
        
        self.glabel = Label(tol=20,tol_type="ppm")
        
        label_tol = wx.StaticText(self.panel, -1, "Tolerance:", pos=(40, base_height))
        
        self.tol_type = "ppm"
        self.tolText = wx.TextCtrl(self.panel, -1, pos=(120,base_height),size=(80,-1),style=wx.ALIGN_RIGHT)
        self.tolText.SetValue("20.0")
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
        
        label_mgf = wx.StaticText(self.panel, -1, 'MGF:', pos=(40,base_height+50))
        self.mgf_file = wx.TextCtrl(self.panel, -1, pos=(120,base_height+50),size=(300,-1))
        self.mgf_file.SetEditable(False)
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
        self.max_mz = wx.TextCtrl(self.panel, -1, pos=(120,base_height+150),size=(80,-1), style=wx.ALIGN_RIGHT, value="4100.0")
        self.activationComboBox = wx.ComboBox(self.panel, -1, value="HCD",
                    pos=(220,base_height+150), choices=["HCD","ETD","ETHCD"],
                    style=wx.CB_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.OnActivationChoose, self.activationComboBox)
        
        self.ShowMassCheckBox = wx.CheckBox(self.panel, -1, label="show mass",
                    pos=(320,base_height+150+5))
        self.ShowMassCheckBox.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.OnMassCheck, self.ShowMassCheckBox)
        
        self.UsepGlycoSite = wx.CheckBox(self.panel, -1, label="pGlycoSite",
                    pos=(420,base_height+150+5))
        self.UsepGlycoSite.SetValue(False)
        self.Bind(wx.EVT_CHECKBOX, self.OnpGlycoSite, self.UsepGlycoSite)
        
        label_spec = wx.StaticText(self.panel, -1, 'Spectrum:', pos=(40,base_height+200))
        self.spec_name = wx.TextCtrl(self.panel, -1, pos=(120,base_height+200),size=(300,-1))
        self.spec_button = wx.Button(self.panel, -1, 'show', pos=(440,base_height+200))
        self.Bind(wx.EVT_BUTTON, self.OnSpecButton, self.spec_button)
        
        box = wx.StaticBox(self.panel, -1, 'self defined glycopeptide', pos=(30, base_height+240), size = (510,120))
        label_glycan = wx.StaticText(self.panel, -1, 'Glycan:', pos=(40, base_height+270))
        self.glycan = wx.TextCtrl(self.panel, -1, pos=(120, base_height+270), size=(300,-1))
        self.glycan.SetValue("H N|8 2|5 2;6 2;7 2;8 2")
        label_peptide = wx.StaticText(self.panel, -1, 'Peptide:', pos=(40, base_height+320))
        self.peptide = wx.TextCtrl(self.panel, -1, pos=(120, base_height+320), size=(300,-1))
        self.peptide.SetValue("ACMJESGK|4|2,Carbamidomethyl[C];3,Oxidation[M]")
        self.pep_button = wx.Button(self.panel, -1, 'show this', pos=(440,base_height+320))
        self.Bind(wx.EVT_BUTTON, self.OnPepButton, self.pep_button)
        
        
        label_batch = wx.StaticText(self.panel, -1, 'Out Folder:', pos=(40,base_height+380))
        self.batch_folder = wx.TextCtrl(self.panel, -1, pos=(120,base_height+380),size=(300,-1))
        self.batch_plot = wx.Button(self.panel, -1, 'Batch Plot', pos=(440,base_height+380))
        self.Bind(wx.EVT_BUTTON, self.OnBatchPlotButton, self.batch_plot)
        
        # self.batch_browse = wx.Button(self.panel, -1, 'Browse:', pos=(40,base_height+380),size=(60,-1))
        # self.Bind(wx.EVT_BUTTON, self.OnBatchBrowseButton, self.batch_browse)
        
    def OnMassCheck(self, event):
        self.glabel.show_mass = self.ShowMassCheckBox.GetValue()
        
    def OnpGlycoSite(self, event):
        self.glabel.use_pGlycoSite = self.UsepGlycoSite.GetValue()
        
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
        openDir = wx.DirDialog(self.panel, "Choose Output Folder","",
                                wx.DD_DEFAULT_STYLE)
        if openDir.ShowModal() == wx.ID_CANCEL:
            return
        self.batch_folder.SetValue(openDir.GetPath())
        
    def OnBatchPlotButton(self, event):
        self.glabel.max_plot_mz = self.GetFloat(self.max_mz.GetValue())
        self.glabel.save_plot(self.batch_folder.GetValue())
        print("****  Finish batch output ****")
    
    def OnMGFButton(self, event):
        openFile=wx.FileDialog(self.panel, "Open MGF file", "", "",
                       "MGF file (*.mgf)|*.mgf|Raw file (*.raw)|*.raw", 
                       wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        if openFile.ShowModal() == wx.ID_CANCEL:
            return
        self.mgf_file.SetValue(openFile.GetPath())
        self.glabel.ReadMGF(self.mgf_file.GetValue())
    
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
        # if config.ResultHasFragments:
        self.glabel.ReadDenovoRes(self.res_file.GetValue())
        
        self.batch_folder.SetValue(os.path.join(os.path.split(self.res_file.GetValue())[0], "gLabel"))
        # else:
            # self.glabel.ReadDBSearchRes(self.res_file.GetValue())
#         if self.type_choose == 1:
#             self.glabel.ReadDenovoRes(self.res_file.GetValue())
#         else:
#             self.glabel.ReadDBSearchRes(self.res_file.GetValue())

class MyApp(wx.App):
    def OnInit(self):
        self.GUI = GUIgLabel(None, wx.ID_ANY, "")
        self.SetTopWindow(self.GUI)
        self.GUI.Show()
        return True
        
        
if __name__ == "__main__":
    print(" *** gLabel in pGlyco ***")
    # ret = os.system("pGlycoDB.exe check_license")
    if True: 
        app = MyApp(0)
        app.MainLoop()
        