import struct
import numpy as np
try:
    import RawFileReader
except Exception:
    RawFileReader = None

ms_reader = {}


def RegisterReader(format, reader):
    ms_reader[format.lower()] = reader


def GetMSReader(filename):
    reader = ms_reader[filename[filename.rfind(".") + 1 :].lower()]()
    reader.open_and_index_file(filename)
    return reader


def get_raw_name(filename):
    return filename[: filename.rfind(".")]


if RawFileReader is not None:
    class RawReader:
        def __init__(self):
            self.file = None
            self.last_scan = 0

        def open_and_index_file(self, filename):
            self.file = RawFileReader(filename)

            self.raw_name = get_raw_name(filename)

            self.last_scan = self.file.GetLastSpectrumNumber()

            self.scanidx = {}
            for scan in range(self.file.GetFirstSpectrumNumber(), self.last_scan + 1):
                if self.file.GetMSOrderForScanNum(scan) == 1:
                    self.scanidx[scan] = (scan, self.file.RTInSecondsFromScanNum(scan))

        def read_a_peaklist(self, scan):
            return self.file.GetCentroidMassListFromScanNum(scan).T

        def close(self):
            if self.file is not None:
                self.file.Close()
                self.file = None

        def __del__(self):
            self.close()

    RegisterReader("raw", RawReader)


class PF1Reader:
    def __init__(self):
        self.file = None
        self.last_scan = 0

    def open_and_index_file(self, filename):
        self.file = open(filename, "rb")

        self.raw_name = get_raw_name(filename)

        self.scanidx = {}
        nspec, ntitle = struct.unpack("2i", self.file.read(8))
        self.file.read(ntitle)
        for i in range(nspec):
            index = self.file.tell()
            scan, nPeak = struct.unpack("2i", self.file.read(8))
            self.file.read(nPeak * 2 * 8 + 4)  # peaks + fake nMix=int
            (RT,) = struct.unpack("d", self.file.read(8))
            self.file.read(4)  # fake charge=int
            self.scanidx[scan] = (index, RT)
            if scan > self.last_scan:
                self.last_scan = scan

    # return (peaklist, RT)
    def read_a_peaklist(self, scan):
        # only read the (mz, inten) list
        self.file.seek(self.scanidx[scan][0])
        scan, nPeak = struct.unpack("2i", self.file.read(8))
        mz_int = struct.unpack(str(nPeak * 2) + "d", self.file.read(nPeak * 2 * 8))
        # n, = struct.unpack("i", self.file.read(4))
        # RT, = struct.unpack("d", self.file.read(8))
        # charge, = struct.unpack("i", self.file.read(4))
        # end scan
        peaklist = []
        for i in range(nPeak):
            mz = mz_int[i * 2]
            inten = mz_int[i * 2 + 1]
            peaklist.append((mz, inten))
        return np.array(peaklist)

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None

    def __del__(self):
        self.close()


RegisterReader("pf1", PF1Reader)


class MS1Reader:
    def __init__(self):
        self.file = None
        self.last_scan = 0

    def open_and_index_file(self, filename):
        self.file = open(filename)

        self.raw_name = get_raw_name(filename)

        self.scanidx = {}
        while True:
            line = self.file.readline()
            if line == "":
                break
            elif line.startswith("S"):
                idx = self.file.tell()
                scan = int(line.split("\t")[1])
                if scan > self.last_scan:
                    self.last_scan = scan
            elif line.startswith("I\tRetTime"):
                RT = float(line.strip().split("\t")[-1])
                self.scanidx[scan] = (idx, RT)

    def read_a_peaklist(self, scan):
        # only read the (mz, inten) list
        self.file.seek(self.scanidx[scan][0])
        peaklist = []
        while True:
            line = self.file.readline()
            if line == "":
                break
            elif line.startswith("S"):
                break
            elif line[0].isdigit():
                mz, inten = line.strip().split()
                mz, inten = float(mz), float(inten)
                peaklist.append((mz, inten))
        return np.array(peaklist)

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None

    def __del__(self):
        self.close()


RegisterReader("ms1", MS1Reader)
