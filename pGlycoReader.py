class GPSM:
    def __init__(self):
        self.spectrum = ""
        self.scan = -1
        self.RT = -1
        self.seq = ""
        self.mod = ""
        self.charge = 0
        self.precursor_mz = 0
        self.glycan = ""
        self.glycan_fragment = []
        self.Y1_glycan = (0,)
        self.glysite = -1
        self.pep_decoy = 0
        self.gly_decoy = 0

        self.pep_score = 0
        self.gly_score = 0
        self.total_score = 0

        self.pep_svmscore = 0
        self.gly_svmscore = 0
        self.total_svmscore = 0

        self.gly_FDR = 0
        self.pep_FDR = 0
        self.total_FDR = 0

        self.pep_feature = []
        self.gly_feature = []

        self.peak_area = -1
        self.mono_area = -1

        self.items = []


def parse_scan_charge_from_spec(spec):
    items = spec.split(".")
    return int(items[-4]), int(items[-2])


class pGlycoReader:
    def __init__(self):
        self.glycan_colname = "Glycan()"
        self.glyco_names = []

    def _set_glyco_name(self):
        self.glyco_names = self.glycan_colname[
            self.glycan_colname.find("(") + 1 : -1
        ].split(",")
        if self.glyco_names[0] == "":
            print(
                "[ERROR] glycan title %s in result file is not correct"
                % self.glycan_colname
            )
            import sys

            sys.exit(-1)

    def _GetGPSM(self, headidx, items):
        psm = GPSM()
        psm.spectrum = items[headidx["GlySpec"]]
        if "Scan" in headidx:
            psm.scan = int(items[headidx["Scan"]])
        else:
            psm.scan, psm.charge = parse_scan_charge_from_spec(self.spectrum)
        psm.seq = items[headidx["Peptide"]]
        psm.mod = items[headidx["Mod"]].strip('"')
        if psm.mod == "null":
            psm.mod = ""
        if "Charge" in headidx:
            psm.charge = int(items[headidx["Charge"]])
        if "PrecursorMZ" in headidx:
            psm.precursor_mz = float(items[headidx["PrecursorMZ"]])
        psm.glycan = items[headidx[self.glycan_colname]].strip(" ").split(" ")
        psm.glycan = [int(i) for i in psm.glycan]
        psm.glysite = int(headidx["GlySite"])
        glyfrag = items[headidx["GlyFrag"]].strip(";")
        if glyfrag:
            glyfrag = glyfrag.split(";")
            for gly in glyfrag:
                psm.glycan_fragment.append(
                    tuple([int(i) for i in gly.strip().split(" ")])
                )
        psm.pep_decoy = int(items[headidx["PepDecoy"]])
        psm.gly_decoy = int(items[headidx["GlyDecoy"]])
        psm.pep_score = float(items[headidx["PepScore"]])
        psm.gly_score = float(items[headidx["GlyScore"]])
        psm.total_score = float(items[headidx["TotalScore"]])
        psm.pep_svmscore = psm.pep_score
        psm.gly_svmscore = psm.gly_score
        psm.total_svmscore = psm.total_score

        if "TotalFDR" in headidx:
            psm.total_FDR = float(items[headidx["TotalFDR"]])
        if "GlycanFDR" in headidx:
            psm.gly_FDR = float(items[headidx["GlycanFDR"]])
        if "PeptideFDR" in headidx:
            psm.pep_FDR = float(items[headidx["PeptideFDR"]])
        if "RT" in headidx:
            psm.RT = float(items[headidx["RT"]])

        return psm

    def ReadAllGPSMs(self, pglyco_file):
        self.gpsm_list = []
        with open(pglyco_file) as f:
            self.head = f.readline().strip().split("\t")
            for h in self.head:
                if h.startswith("Glycan("):
                    self.glycan_colname = h
                    self._set_glyco_name()
                    break
            self.headidx = dict(zip(self.head, range(len(self.head))))
            lines = f.readlines()
            for line in lines:
                items = line.strip().split("\t")
                if items[self.headidx["Rank"]] != "1":
                    continue
                self.gpsm_list.append(self._GetGPSM(self.headidx, items))
                self.gpsm_list[-1].items = items
        return self.gpsm_list

    def WriteAllGPSMs(self, pglyco_file):
        with open(pglyco_file, "w") as f:
            f.write("\t".join(self.head) + "\n")
            for gpsm in self.gpsm_list:
                f.write("\t".join(gpsm.items))
                f.write("\n")
