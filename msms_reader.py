import os
import struct
import numpy as np
from RawFileReader import RawFileReader

def get_scan_charge_by_specname(specname): #pParse format
    items = specname.split(".")
    return int(items[-4]), int(items[-3]) #scan, charge
    
def get_raw_name(filename):
    filename = os.path.basename(filename)
    if filename.rfind("_") != -1:
        return filename[:filename.rfind("_")]
    else:
        return filename[:filename.rfind(".")]
    
msms_reader = {}
def RegisterReader(format, reader):
    msms_reader[format.lower()] = reader
        
def GetMSMSReader(filename):
    reader = msms_reader[filename[filename.rfind(".")+1:].lower()]()
    reader.open_and_index_file(filename)
    return reader
    
class RawReader:
    def __init__(self):
        self.file = None
        self.last_scan = 0
    
    def open_and_index_file(self, filename):
        self.file = RawFileReader(filename)
        
        self.raw_name = get_raw_name(filename)
        
        self.last_scan = self.file.GetLastSpectrumNumber()
        
        self.scanidx = {}
        for scan in range(self.file.GetFirstSpectrumNumber(), self.last_scan+1):
            if self.file.GetMSOrderForScanNum(scan) == 2:
                self.scanidx[scan] = (scan, self.file.RTInSecondsFromScanNum(scan))
        
    def read_a_peaklist(self, scan):
        return self.file.GetCentroidMassListFromScanNum(scan).T
        
    def close(self):
        if self.file is not None:
            self.file.Close()
            self.file = None
            
    def __del__(self):
        self.close()
RegisterReader('raw', RawReader)

class PF2Reader:
    def __init__(self):
        self.file = None
        
    def open_and_index_file(self, filename):
        self.file = open(filename, 'rb')
        
        self.raw_name = get_raw_name(filename)
        
        self.scanidx = {}
        if False and os.path.isfile(filename+'idx'):
            f = open(filename+'idx','rb')
            while True:
                chunk = f.read(8)
                if not chunk: break
                scan, index = struct.unpack('2i',chunk)
                print(scan)
                self.scanidx[scan] = index
            f.close()
        else:
            nspec,ntitle = struct.unpack("2i",self.file.read(8))
            self.file.read(ntitle)
            for i in range(nspec):
                index = self.file.tell()
                scan, nPeak = struct.unpack("2i",self.file.read(8))
                self.file.read(nPeak*2*8) #peaks
                nMix, = struct.unpack("i", self.file.read(4))
                self.file.read(nMix*12)
                self.scanidx[scan] = index
                
    def read_a_peaklist(self, scan):
        # only read the (mz, inten) list
        self.file.seek(self.scanidx[scan])
        scan, nPeak = struct.unpack("2i",self.file.read(8))
        mz_int = struct.unpack(str(nPeak*2)+"d", self.file.read(nPeak*2*8))
        peaklist = []
        for i in range(nPeak):
            mz = mz_int[i*2]
            inten = mz_int[i*2+1]
            peaklist.append( (mz, inten) )
        return np.array(peaklist)
        
    def close(self):
        if self.file is not None: 
            self.file.close()
            self.file = None
        
    def __del__(self):
        self.close()
RegisterReader("pf2", PF2Reader)
        
class MS2Reader:
    def __init__(self):
        self.file = None
        
    def open_and_index_file(self, filename):
        self.file = open(filename)
        
        self.raw_name = get_raw_name(filename)
        
        self.scanidx = {}
        while True:
            line = self.file.readline()
            if line == "": break
            elif line.startswith("S"):
                idx = self.file.tell()
                scan = int(line.split("\t")[1])
                self.scanidx[scan] = idx
    
    def read_a_peaklist(self, scan):
        # only read the (mz, inten) list
        self.file.seek(self.scanidx[scan])
        peaklist = []
        while True:
            line = self.file.readline()
            if line == "": break
            elif line.startswith("S"): break
            elif line[0].isdigit():
                mz, inten = line.strip().split()
                mz, inten = float(mz), float(inten)
                peaklist.append( (mz, inten) )
        return np.array(peaklist)
        
    def close(self):
        if self.file is not None: 
            self.file.close()
            self.file = None
        
    def __del__(self):
        self.close()
RegisterReader("ms2", MS2Reader)
        
class MGFReader:
    def __init__(self, pParse = True):
        self.file = None
        self.pParse_format = pParse
        
    def open_and_index_file(self, filename):
        self.file = open(filename)
        
        self.raw_name = get_raw_name(filename)
        
        self.scanidx = {}
        while True:
            line = self.file.readline()
            if line == "": break
            elif line.startswith("BEGIN IONS"):
                scan = -1
                idx = self.file.tell()
            elif line.startswith("SCAN"):
                if scan != -1:
                    scan = int(line[line.find("=")+1:])
                    self.scanidx[scan] = idx
            elif line.startswith("TITLE"):
                if self.pParse_format:
                    scan,_ = get_scan_charge_by_specname(line.strip())
                    self.scanidx[scan] = idx
    
    def read_a_peaklist(self, scan):
        # only read the (mz, inten) list
        self.file.seek(self.scanidx[scan])
        peaklist = []
        while True:
            line = self.file.readline()
            if line == "": break
            elif line.startswith("END IONS"): break
            elif line[0].isdigit():
                mz, inten = line.strip().split()
                mz, inten = float(mz), float(inten)
                peaklist.append( (mz, inten) )
        return np.array(peaklist)
        
    def close(self):
        if self.file is not None: 
            self.file.close()
            self.file = None
        
    def __del__(self):
        self.close()
RegisterReader("mgf", MGFReader)
