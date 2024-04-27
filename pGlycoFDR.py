from FDR_estimator import FMM_FDR, TDA_FDR
from pGlyco_config import pGlycoConfig

import os
import sys

short_glycan_len = 4
verbose = False
FMM_verbose = verbose


def PPM2Da(mz, ppm):
    return ppm * mz / 1e6


class pGlycoFilterMix:
    def __init__(self):
        self.isotope_ppm = 20.0

        self.pParse_unreliable_mix = [
            # ([high_mass_list], [low_mass_list], delta_mass)
            ([("F", 2)], [("A", 1)], 1.003),  # 2F = 1A + 1.003, remove 2*F
            ([("F", 4)], [("A", 2)], 2.006),  # 4F = 2A + 2.006, remove 2*F
            (
                [("H", 2), ("N", 2)],
                [("F", 3), ("A", 1)],
                1.003,
            ),  # 2H + 2N = 3F + 1A + 1.003, remove 2H + 2N
        ]

    def Filter_unreliable_pParse(self, pglyco_table, score_factor=10000):
        orig_res = {}
        for res in pglyco_table._data:
            # scan = ExtractScanFromSpec(self.GetItemByColName(res, "PepSpec"))
            scan = pglyco_table.GetItemByColName(res, "PepSpec")
            scan = ".".join(scan.split(".")[:-3])
            if scan in orig_res:
                orig_res[scan].append(res)
            else:
                orig_res[scan] = [res]
        pglyco_table._data = []

        def delta_glycan(onescan, i, j):
            glycani = (
                pglyco_table.GetItemByColName(onescan[i], pglyco_table.GlycanCol)
                .strip()
                .split(" ")
            )
            glycanj = (
                pglyco_table.GetItemByColName(onescan[j], pglyco_table.GlycanCol)
                .strip()
                .split(" ")
            )
            glycani = [int(glyco) for glyco in glycani]
            glycanj = [int(glyco) for glyco in glycanj]
            return [glycani[t] - glycanj[t] for t in range(len(pglyco_table.GlycoIdx))]

        def check_can_remove(glycan_delta, high_mass_gly, low_mass_gly):
            for x in high_mass_gly:
                if glycan_delta[pglyco_table.GlycoIdx[x[0]]] != x[1]:
                    return False
            for x in low_mass_gly:
                if glycan_delta[pglyco_table.GlycoIdx[x[0]]] != -x[1]:
                    return False
            return True

        for onescan in orig_res.values():
            removeTag = [0] * len(onescan)
            for i in range(len(onescan)):
                for j in range(i + 1, len(onescan)):
                    if pglyco_table.GetItemByColName(
                        onescan[i], "Peptide"
                    ) != pglyco_table.GetItemByColName(onescan[j], "Peptide"):
                        continue
                    if pglyco_table.GetItemByColName(
                        onescan[i], "Charge"
                    ) != pglyco_table.GetItemByColName(onescan[j], "Charge"):
                        continue

                    MHi = float(
                        pglyco_table.GetItemByColName(onescan[i], "PrecursorMH")
                    )
                    MHj = float(
                        pglyco_table.GetItemByColName(onescan[j], "PrecursorMH")
                    )

                    PepScorei = float(
                        pglyco_table.GetItemByColName(onescan[i], "PepScore")
                    )
                    PepScorej = float(
                        pglyco_table.GetItemByColName(onescan[j], "PepScore")
                    )

                    GlyScorei = float(
                        pglyco_table.GetItemByColName(onescan[i], "GlyScore")
                    )
                    GlyScorej = float(
                        pglyco_table.GetItemByColName(onescan[j], "GlyScore")
                    )

                    PPMi = abs(float(pglyco_table.GetItemByColName(onescan[i], "PPM")))
                    PPMj = abs(float(pglyco_table.GetItemByColName(onescan[j], "PPM")))

                    if pglyco_table.GetItemByColName(
                        onescan[i], "Mod"
                    ) != pglyco_table.GetItemByColName(onescan[j], "Mod"):
                        if PepScorei > PepScorej:
                            removeTag[j] = 1
                        elif PepScorei < PepScorej:
                            removeTag[i] = 1
                        else:
                            if MHi > MHj:
                                removeTag[i] = 1
                            else:
                                removeTag[j] = 1
                    else:
                        if (
                            abs(abs(MHi - MHj) - round(abs(MHi - MHj))) > 0.05
                        ):
                            if GlyScorei < GlyScorej:
                                removeTag[i] = 1
                            else:
                                removeTag[j] = 1
                        elif MHi > MHj:
                            removeTag[i] = 1
                        else:
                            removeTag[j] = 1

                    # for high_mass_gly, low_mass_gly, delta_mass in self.pParse_unreliable_mix:
                    # no_high_low = False
                    # for x in high_mass_gly:
                    # if x[0] not in pglyco_table.GlycoIdx: no_high_low = True
                    # for x in low_mass_gly:
                    # if x[0] not in pglyco_table.GlycoIdx: no_high_low = True
                    # if no_high_low: continue

                    ## if MHi is higher
                    # if abs(MHi - MHj - delta_mass) <= PPM2Da(MHi, self.isotope_ppm):
                    # glycan_delta = delta_glycan(onescan, i, j)
                    # if check_can_remove(glycan_delta, high_mass_gly, low_mass_gly) and float(pglyco_table.GetItemByColName(onescan[i], "GlyScore")) < float(pglyco_table.GetItemByColName(onescan[j], "GlyScore"))*score_factor:
                    # removeTag[i] = 1
                    # break

                    ## if MHj is higher
                    # elif abs(MHj - MHi - delta_mass) <= PPM2Da(MHj, self.isotope_ppm):
                    # glycan_delta = delta_glycan(onescan, j, i)
                    # if check_can_remove(glycan_delta, high_mass_gly, low_mass_gly) and float(pglyco_table.GetItemByColName(onescan[j], "GlyScore")) < float(pglyco_table.GetItemByColName(onescan[i], "GlyScore"))*score_factor:
                    # removeTag[j] = 1
                    # break

            for i in range(len(removeTag)):
                if removeTag[i] == 0:
                    pglyco_table._data.append(onescan[i])


class pGlycoFDR:
    def __init__(self):
        self.pep_use_FMM = 0

    def FDREstimate(self, pglyco_table):
        if len(pglyco_table.targetPepScore) < 10:
            pglyco_table.AddValues_IsShortGlycan()
            pglyco_table.AddValues([0] * len(pglyco_table.targetPepScore), "GlycanPEP")
            pglyco_table.AddValues([0] * len(pglyco_table.targetPepScore), "GlycanFDR")
            pglyco_table.AddValues([0] * len(pglyco_table.targetPepScore), "PeptidePEP")
            pglyco_table.AddValues([0] * len(pglyco_table.targetPepScore), "PeptideFDR")
            pglyco_table.AddValues([0] * len(pglyco_table.targetPepScore), "TotalFDR")
            return
        if "TotalFDR" in pglyco_table.headeridx:
            return

        # peptide
        if self.pep_use_FMM:
            self.PepFDR = FMM_FDR()
        else:
            self.PepFDR = TDA_FDR()
        self.PepFDR.fit(
            pglyco_table.targetPepScore, pglyco_table.decoyPepScore, verbose=FMM_verbose
        )
        pepfdr, peppep = self.PepFDR.estimate(pglyco_table.GetColumn("PepScore", float))

        # decoy
        self.GlyFDR = FMM_FDR()
        self.GlyFDR.fit(
            pglyco_table.targetGlyScore, pglyco_table.decoyGlyScore, verbose=FMM_verbose
        )
        glyscores = pglyco_table.GetColumn("GlyScore", float)

        # G cap P
        self.CapFDR = FMM_FDR()
        self.CapFDR.fit(
            pglyco_table.targetTotalScore,
            pglyco_table.decoyTotalScore,
            verbose=FMM_verbose,
        )
        totalscores = pglyco_table.GetColumn("TotalScore", float)

        not_short_glycan_idxes = []
        _glyscores = []
        _totalscores = []
        for i in range(len(totalscores)):
            if not pglyco_table.IsShortGlycan(pglyco_table._data[i]):
                not_short_glycan_idxes.append(i)
                _totalscores.append(totalscores[i])
                _glyscores.append(glyscores[i])

        _glyfdr, _glypep = self.GlyFDR.estimate(_glyscores)
        _totalfdr, totalpep = self.CapFDR.estimate(_totalscores)
        glyfdr = [0] * len(totalscores)
        glypep = [-1] * len(totalscores)
        totalfdr = [0] * len(totalscores)
        for i in range(len(not_short_glycan_idxes)):
            glyfdr[not_short_glycan_idxes[i]] = _glyfdr[i]
            glypep[not_short_glycan_idxes[i]] = _glypep[i]
            totalfdr[not_short_glycan_idxes[i]] = _totalfdr[i]

        for i in range(len(totalfdr)):
            fdr = glyfdr[i] + pepfdr[i] - totalfdr[i]
            totalfdr[i] = min(max((glyfdr[i], pepfdr[i], fdr)), 1.0)
        pglyco_table.AddValues_IsShortGlycan()
        pglyco_table.AddValues(glypep, "GlycanPEP")
        pglyco_table.AddValues(glyfdr, "GlycanFDR")
        pglyco_table.AddValues(peppep, "PeptidePEP")
        pglyco_table.AddValues(pepfdr, "PeptideFDR")
        pglyco_table.AddValues(totalfdr, "TotalFDR")

    def FMM_test(self, table):
        import matplotlib.pyplot as plt

        self.GlyFDR.get_fmm_model().plot(
            "glycan", table.GetColumn("GlyScore", float), table.decoyGlyScore
        )
        if self.pep_use_FMM:
            self.PepFDR.get_fmm_model().plot(
                "peptide", table.GetColumn("PepScore", float), table.decoyPepScore
            )
        self.TotalFDR.get_fmm_model().plot(
            "g-cap-p", table.GetColumn("TotalScore", float), table.decoyTotalScore
        )

        plt.show()


class pGlycoTable:
    def __init__(self):
        self.reset()

    def reset(self):
        self.decoyGlyScore = []
        self.decoyPepScore = []
        self.decoyPepDict = {}
        self.decoyTotalScore = []

        self.targetPepScore = []
        self.targetGlyScore = []
        self.targetTotalScore = []

        self.headeridx = {}
        self.header = ""

        self._data = []
        self._pep_decoy_data = []
        self._gly_decoy_data = []
        self._all_decoy_data = []

        self.GlycanCol = "Glycan"
        self.GlycoIdx = dict(zip("HNAGF", range(len("HNAGF"))))

    # cannot estimate FDR for short glycans
    def IsShortGlycan(self, row_items):
        if "PlausibleStruct" in self.headeridx:
            return (
                self.GetItemByColName(row_items, "PlausibleStruct").count(")")
                <= short_glycan_len
            )
        else:
            return (
                sum(
                    [
                        int(i)
                        for i in self.GetItemByColName(row_items, self.GlycanCol)
                        .strip()
                        .split(" ")
                    ]
                )
                <= short_glycan_len
            )

    def RearangeGlycoIdx(self):
        glycos = self.GlycanCol[
            self.GlycanCol.find("(") + 1 : self.GlycanCol.find(")")
        ].split(",")
        self.GlycoIdx = dict(zip(glycos, range(len(glycos))))

    def AddValues_IsShortGlycan(self):
        if "IsSmallGlycan" not in self.headeridx:
            self.AddValues(
                ["1" if self.IsShortGlycan(items) else "0" for items in self._data],
                "IsSmallGlycan",
            )

    def AddColumn(self, colname):
        self.headeridx[colname] = len(self.header)
        self.header.append(colname)

    def AddValues(self, values, colname):
        self.AddColumn(colname)
        for i in range(len(values)):
            self._data[i].append(str(values[i]))

    def GetColumn(self, colname, dtype=float, nonredundant=None):
        if nonredundant is None:
            return [dtype(self.GetItemByColName(data, colname)) for data in self._data]
        else:
            data_dict = {}
            for data in self._data:
                key = self.GetItemByColName(data, nonredundant)
                val = dtype(self.GetItemByColName(data, colname))
                if key in data_dict:
                    if data_dict[key] < val:
                        data_dict[key] = val
                else:
                    data_dict[key] = val
            return data_dict

    def GetItemByColName(self, items, colname):
        if colname in self.headeridx:
            return items[self.headeridx[colname]]
        else:
            return None

    def WriteTable(self, filename, FDR):
        print(filename)
        f = open(filename, "w")
        f.write("\t".join(self.header) + "\n")
        fdr_count = 0
        for item in self._data:
            if float(self.GetItemByColName(item, "TotalFDR")) <= FDR:
                f.write("\t".join(item) + "\n")
                fdr_count += 1
        f.close()
        print("[pGlycoFDR] %d GPSMs at %.1f%% FDR\n" % (fdr_count, FDR * 100))

    def ReadTable(self, filename):
        self.reset()
        try:
            f = open(filename)
        except:
            print("Could not open %s!" % filename)
            sys.exit(-1)
        line = f.readline()
        items = line.strip().split("\t")
        for item in items:
            if item.startswith("Glycan("):
                self.GlycanCol = item
                self.RearangeGlycoIdx()
                break
        self.header = items
        self.headeridx = dict(zip(items, range(len(items))))
        while True:
            line = f.readline()
            if line == "":
                break
            items = line.strip().split("\t")
            if self.GetItemByColName(items, "Rank") == "1":
                if (
                    self.GetItemByColName(items, "GlyDecoy") == "1"
                    and self.GetItemByColName(items, "PepDecoy") == "1"
                ):
                    if not self.IsShortGlycan(items):
                        self.decoyTotalScore.append(
                            float(self.GetItemByColName(items, "TotalScore"))
                        )
                    self._all_decoy_data.append(items)
                elif self.GetItemByColName(items, "GlyDecoy") == "1":
                    if not self.IsShortGlycan(items):
                        self.decoyGlyScore.append(
                            float(self.GetItemByColName(items, "GlyScore"))
                        )
                    self._gly_decoy_data.append(items)
                elif self.GetItemByColName(items, "PepDecoy") == "1":
                    score = float(self.GetItemByColName(items, "PepScore"))
                    self.decoyPepScore.append(score)

                    # modseq = self.GetItemByColName(items, "Peptide")+self.GetItemByColName(items, "Mod")
                    # if modseq in self.decoyPepDict:
                    # if score > self.decoyPepDict[modseq]: self.decoyPepDict[modseq] = score
                    # else:
                    # self.decoyPepDict[modseq] = score

                    # self._pep_decoy_data.append(items)
                else:
                    if not self.IsShortGlycan(items):
                        self.targetGlyScore.append(
                            float(self.GetItemByColName(items, "GlyScore"))
                        )
                        self.targetTotalScore.append(
                            float(self.GetItemByColName(items, "TotalScore"))
                        )
                    self.targetPepScore.append(
                        float(self.GetItemByColName(items, "PepScore"))
                    )
                    self._data.append(items)
        f.close()
        if len(self.decoyPepDict) > 0:
            self.decoyPepScore = list(self.decoyPepDict.values())
        if verbose:
            print("[pGlycoFDR] %d PSMs are from both target" % len(self._data))
        if verbose:
            print("[pGlycoFDR] %d PSMs are from glycan decoy" % len(self.decoyGlyScore))
        if verbose:
            print(
                "[pGlycoFDR] %d PSMs are from peptide decoy " % len(self.decoyPepScore)
            )
        if verbose:
            print("[pGlycoFDR] %d PSMs are from both decoy" % len(self.decoyTotalScore))


def Run(filename, config, test="0"):
    table = pGlycoTable()
    table.ReadTable(filename)
    pGlycoFilterMix().Filter_unreliable_pParse(table)
    m = pGlycoFDR()
    m.pep_use_FMM = config.FMM_for_peptide_FDR
    m.FDREstimate(table)
    if test != "0":
        m.FMM_test(table)
    table.WriteTable(filename[:-4] + "-FDR.txt", config.FDR)
    table.WriteTable(filename[:-4] + "-FDR-noFiltered.txt", 1)


def ParseResults(config):
    ret = []
    if config.percolator and os.path.isfile(
        os.path.join(config.output_dir, "%s-GP-Raw1-SVM.txt" % config.pGlycoType)
    ):
        suff = "-SVM"
    else:
        suff = ""
    for i in range(len(config.mgf_list)):
        ret.append(
            os.path.join(
                config.output_dir,
                "%s-GP-Raw%d%s.txt" % (config.pGlycoType, i + 1, suff),
            )
        )
    return ret


def MergeFDRFile(fdr_files, pGlycoType):
    f = open(fdr_files[0])
    lines = f.readlines()
    f.close()
    outfolder = os.path.split(fdr_files[0])[0]
    outfile = os.path.join(outfolder, "%s-GP-FDR.txt" % pGlycoType)

    print("merge into %s" % outfile)
    out = open(outfile, "w")
    out.writelines(lines)
    for i in range(1, len(fdr_files)):
        f = open(fdr_files[i])
        lines = f.readlines()
        f.close()
        out.writelines(lines[1:])
    out.close()


def ParseArgs(argv):
    arg_dict = {}
    i = 1
    while i < len(argv):
        arg_dict[argv[i]] = argv[i + 1]
        i += 2
    return arg_dict


if __name__ == "__main__":
    argd = ParseArgs(sys.argv)

    def get_test():
        if "-test" in argd:
            return argd["-test"]
        else:
            return "0"

    pGlycoType = "pGlycoDB"  # or "pGlycoNovo"
    if "-t" in argd:
        pGlycoType = argd["-t"]
    if "-r" in argd:
        Run(argd["-r"], pGlycoConfig(), get_test())
    elif "-p" in argd:
        config = pGlycoConfig(argd["-p"])
        files = ParseResults(config)
        fdr_files = []
        for filename in files:
            Run(filename, config, get_test())
            fdr_files.append(filename[:-4] + "-FDR.txt")
        MergeFDRFile(fdr_files, config.pGlycoType)
    else:
        print("Error argument!")
        print("Usage: pGlycoFDR option filename")
        print("  -p: filename = parameter file")
        print("  -r: filename = pGlyco result file")
