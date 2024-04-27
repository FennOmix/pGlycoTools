from isotope import IsotopeSimple as isotp
from ms_reader import GetMSReader
from pGlyco_config import pGlycoConfig
import spectrum_index as spi
import quant
from pGlycoReader import pGlycoReader
from mass_utils import GlyIonCalc, PepIonCalc

import struct
import numpy as np

import os


mass_proton = 1.007276
mass_isotope = 1.0033
tol = 20.0
ppm = True

marker_ions = np.array([273.0848518414 + mass_proton, 291.09541652769997 + mass_proton])
mixed_delta = 146.0579088094 * 2 - 291.09541652769997


class pf2reader:
    def __init__(self, pf2=None):
        self.scanidx = {}
        self.pf2 = None
        if pf2:
            self.open(pf2)

    def close(self):
        if self.pf2 is not None:
            self.pf2.close()
            self.pf2 = None

    def open(self, pf2):
        self.close()
        self.pf2 = open(pf2, "rb")

        self.scanidx = {}
        with open(pf2 + "idx", "rb") as f:
            while True:
                chunk = f.read(8)
                if not chunk:
                    break
                scan, index = struct.unpack("2i", chunk)
                self.scanidx[scan] = index

    def read_peaklist(self, scan):
        # only read the (mz, inten) list
        if scan not in self.scanidx:
            return None, None
        self.pf2.seek(self.scanidx[scan])
        _scan, nPeak = struct.unpack("2i", self.pf2.read(8))
        mz_int = struct.unpack(str(nPeak * 2) + "d", self.pf2.read(nPeak * 2 * 8))
        masses = []
        intens = []
        for i in range(nPeak):
            mz = mz_int[i * 2]
            inten = mz_int[i * 2 + 1]
            masses.append(mz)
            intens.append(inten)
        return np.array(masses), np.array(intens)


def match(masses, intens, query_masses):
    query_masses = np.sort(query_masses)
    if ppm:
        tols = query_masses * tol * 1e-6
    else:
        tols = np.ones_like(query_masses) * tol
    i = 0
    j = 0
    query_intens = np.zeros_like(query_masses)
    while i < len(masses) and j < len(query_masses):
        if masses[i] > (query_masses[j] + tols[j]):
            j += 1
        elif masses[i] < (query_masses[j] - tols[j]):
            i += 1
        else:
            query_intens[j] += intens[i]
            i += 1
    return query_intens


class XICExtractor:
    def __init__(self, cfg, smooth="savgol_filter", window=21):
        self.config = pGlycoConfig(cfg)
        self.isppm = self.config.pre_isppm
        self.tol = self.config.pre_tol + 0.003
        self.RT_win = self.config.MS1_RT_win / 2
        self.Scan_win = int(self.config.MS1_Scan_win / 2)
        self.RT_group = 60
        self.Scan_group = 40

        self.inter_thres = 1.5
        self.glycalc = None
        self.pepcalc = PepIonCalc(self.config.mod_dict, self.config.aa_dict)
        self.isotope_calc = isotp(self.config)

        self.quant = quant.Quant_smooth(window)

    def _get_ms_reader(self, rawname):
        if os.path.isfile(rawname + ".pf1"):
            print("[XIC] Indexing %s.pf1" % rawname)
            return GetMSReader(rawname + ".pf1")
        elif os.path.isfile(rawname + ".ms1"):
            print("[XIC] Indexing %s.ms1" % rawname)
            return GetMSReader(rawname + ".ms1")
        else:
            print('[XIC] ms1 file of raw "%s" not found' % rawname)
            return None

    def _get_pf2(self, mgf):
        if mgf.lower().endswith(".raw"):
            if os.path.isfile(mgf[:-4] + "_HCDFT.pf2"):
                return mgf[:-4] + "_HCDFT.pf2"
            elif os.path.isfile(mgf[:-4] + "_ETHCDFT.pf2"):
                return mgf[:-4] + "_ETHCDFT.pf2"
        elif mgf[:-4].endswith("_HCDpdETXXDFT"):
            return mgf[: -len("_HCDpdETXXDFT.mgf")] + "_HCDFT.pf2"
        elif mgf[:-4].endswith("_ETDFT"):
            return mgf[: -len("_ETDFT.mgf")] + "_ETDFT.pf2"
        elif mgf[:-4].endswith("_ETHCDFT"):
            return mgf[: -len("_ETHCDFT.mgf")] + "_ETHCDFT.pf2"
        elif mgf[:-4].endswith("_HCDFT"):
            return mgf[: -len("_HCDFT.mgf")] + "_HCDFT.pf2"
        else:
            print("Unknow pf2 for mgf '{}'".format(mgf))
            return None

    def _get_pf2_reader(self, mgf):
        pf2 = self._get_pf2(mgf)
        if pf2 is None or not os.path.isfile(pf2):
            return None
        return pf2reader(pf2)

    def _get_rt(self, gpsm):
        return gpsm.RT

    def _de_interference(self, isotope_dist, matched_inten):
        if matched_inten[0] == 0:
            return np.zeros(len(isotope_dist))
        less_inter_idx = 0
        for i in range(1, len(isotope_dist)):
            if isotope_dist[i] <= 0.1 or matched_inten[i] == 0:
                break
            theo_ratio = isotope_dist[i] / isotope_dist[i - 1]
            real_ratio = matched_inten[i] / matched_inten[i - 1]
            if real_ratio > theo_ratio * self.inter_thres:
                break
            elif real_ratio < theo_ratio / self.inter_thres:
                less_inter_idx = i
        return (
            matched_inten[less_inter_idx] / isotope_dist[less_inter_idx] * isotope_dist
        )

    def _check_window(self, rt, scan, gpsm_rt, gpsm_scan, rt_win, scan_win):
        if gpsm_rt > 0:
            return abs(rt - gpsm_rt) <= rt_win
        else:
            return abs(scan - gpsm_scan) <= scan_win

    def _out_win_right(self, rt, scan, gpsm_rt, gpsm_scan, rt_win, scan_win):
        if gpsm_rt > 0:
            return rt > gpsm_rt + rt_win
        else:
            return scan - gpsm_scan <= scan_win

    def group_by_raw_RT(self, gpsm_list, headidx):
        gpsm_group = {}
        sort_list = sorted(gpsm_list, key=lambda x: x.scan)
        for gpsm in sort_list:
            rawname = gpsm.items[headidx["RawName"]].lower()
            if rawname in gpsm_group:
                _prev_gpsm = gpsm_group[rawname][-1][-1]
                if self._check_window(
                    gpsm.RT,
                    gpsm.scan,
                    _prev_gpsm.RT,
                    _prev_gpsm.scan,
                    self.RT_group,
                    self.Scan_group,
                ):
                    gpsm_group[rawname][-1].append(gpsm)
                else:
                    gpsm_group[rawname].append([gpsm])
            else:
                gpsm_group[rawname] = [[gpsm]]
        return gpsm_group

    def _get_rawname(self, mgf):
        idx = mgf[:-4].rfind("_" + self.config.fragmentation)
        if idx == -1:
            if mgf[:-4].endswith("_HCDpdETXXDFT"):
                return mgf[: -len("_HCDpdETXXDFT.mgf")]
            elif mgf[:-4].endswith("_ETDFT"):
                return mgf[: -len("_ETDFT.mgf")]
            elif mgf[:-4].endswith("_ETHCDFT"):
                return mgf[: -len("_ETHCDFT.mgf")]
            elif mgf[:-4].endswith("_HCDFT"):
                return mgf[: -len("_HCDFT.mgf")]
            else:
                return mgf[: mgf.rfind(".")]
        else:
            return mgf[:idx]

    def PeakArea(self, pglyco_file=None):
        pglyco_file = (
            os.path.join(
                self.config.output_dir, self.config.pGlycoType + "-GP-FDR-Pro.txt"
            )
            if pglyco_file is None
            else pglyco_file
        )

        print("[XIC] Loading pGlyco results ...")
        gpsm_reader = pGlycoReader()
        tmp_gpsm_list = gpsm_reader.ReadAllGPSMs(pglyco_file)

        self.config.set_glyco_names(gpsm_reader.glyco_names)
        if "F" in gpsm_reader.glyco_names and "A" in gpsm_reader.glyco_names:
            fuc_idx = gpsm_reader.glyco_names.index("F")
            ac_idx = gpsm_reader.glyco_names.index("A")
        else:
            fuc_idx = None
            ac_idx = None
        # glycan_col_idx = gpsm_reader.head.index(gpsm_reader.glycan_colname)
        # glycan_com_idx = gpsm_reader.head.index("GlycanComposition")
        self.glycalc = GlyIonCalc(gpsm_reader.glyco_names, self.config.gly_dict)
        if "Ion_274.09" in gpsm_reader.headidx:
            ion_274_idx = gpsm_reader.headidx["Ion_274.09"]
        else:
            ion_274_idx = None
        if "Ion_292.10" in gpsm_reader.headidx:
            ion_292_idx = gpsm_reader.headidx["Ion_292.10"]
        else:
            ion_292_idx = None
        if "Ion_204.09" in gpsm_reader.headidx:
            ion_204_idx = gpsm_reader.headidx["Ion_204.09"]
        else:
            ion_204_idx = None

        print("[XIC] RT window is [-%.1f, +%.1f] seconds" % (self.RT_win, self.RT_win))

        reader_dict = {}
        # pf2reader_dict = {}
        for mgf_file in self.config.mgf_list:
            rawname = self._get_rawname(mgf_file).lower()
            reader_dict[os.path.basename(rawname)] = self._get_ms_reader(rawname)
            # pf2reader_dict[os.path.basename(rawname)] = self._get_pf2_reader(mgf_file)

        # fqd = open(pglyco_file[:-4]+"-Quant-Data.txt","w",32*1024*1024)
        fq = open(pglyco_file[:-4] + "-Quant.txt", "w", 32 * 1024 * 1024)
        # gpsm_reader.head.append("LabeledMz")
        gpsm_reader.head.append("MonoArea")
        gpsm_reader.head.append("IsotopeArea")
        gpsm_reader.head.append("Corrected" + gpsm_reader.glycan_colname)
        gpsm_reader.head.append("CorrectedComposition")
        gpsm_reader.head.append("CorrectedMonoArea")
        gpsm_reader.head.append("CorrectedIsotopeArea")

        empty_correct = ["", "", "-1", "-1"]
        empty_head = ["-1", "-1"] + empty_correct

        fq.write("\t".join(gpsm_reader.head) + "\n")

        gpsm_group = self.group_by_raw_RT(tmp_gpsm_list, gpsm_reader.headidx)
        cnt = 0
        for rawname, one_raw_list in gpsm_group.items():
            if rawname in reader_dict:
                _reader = reader_dict[rawname]
            else:
                _reader = None
            # if rawname in pf2reader_dict: _pf2reader = pf2reader_dict[rawname]
            # else: _pf2reader = None
            if not _reader:
                for gpsm_list in one_raw_list:
                    for gpsm in gpsm_list:
                        fq.write("\t".join(gpsm.items))
                        fq.write("\t%s\n" % ("\t".join(empty_head)))
                    cnt += len(gpsm_list)
                    print(
                        "[XIC] Quantifying: %.1f%% (%d/%d)"
                        % (
                            float(cnt) / len(tmp_gpsm_list) * 100,
                            cnt,
                            len(tmp_gpsm_list),
                        ),
                        end="\r",
                    )
                continue
            # with open(os.path.join(self.config.output_dir, "ScanRT.txt"),"w") as f:
            # for scan, (idx, RT) in _reader.scanidx.items():
            # f.write("%d\t%d\t%f\n"%(scan, idx, RT))
            for gpsm_list in one_raw_list:
                ms1_list = []
                scan = gpsm_list[0].scan - 1
                while scan > 0:
                    if scan in _reader.scanidx:
                        RT = _reader.scanidx[scan][1]
                        if not self._check_window(
                            RT,
                            scan,
                            gpsm_list[0].RT,
                            gpsm_list[0].scan,
                            self.RT_win,
                            self.Scan_win,
                        ):
                            break
                        peaklist = _reader.read_a_peaklist(scan)
                        peakidx = spi.PeakIndexing(peaklist, self.tol, self.isppm)
                        ms1_list.append((peaklist, peakidx, RT, scan))
                    scan -= 1
                ms1_list = ms1_list[::-1]

                scan = gpsm_list[0].scan + 1
                while scan <= _reader.last_scan:
                    if scan in _reader.scanidx:
                        RT = _reader.scanidx[scan][1]
                        if self._out_win_right(
                            RT,
                            scan,
                            gpsm_list[-1].RT,
                            gpsm_list[-1].scan,
                            self.RT_win,
                            self.Scan_win,
                        ):
                            break
                        peaklist = _reader.read_a_peaklist(scan)
                        peakidx = spi.PeakIndexing(peaklist, self.tol, self.isppm)
                        ms1_list.append((peaklist, peakidx, RT, scan))
                    scan += 1

                def search_precursor_scan(gpsm_list, ms1_list):
                    i = 0
                    j = 1
                    precursor_idx = {}
                    while i < len(gpsm_list) and j < len(ms1_list):
                        if ms1_list[j - 1][-1] < gpsm_list[i].scan:
                            if ms1_list[j][-1] >= gpsm_list[i].scan:
                                precursor_idx[gpsm_list[i].scan] = j - 1
                                i += 1
                            else:
                                j += 1
                        else:
                            i += 1
                    for gpsm in gpsm_list:
                        if gpsm.scan not in precursor_idx:
                            precursor_idx[gpsm.scan] = len(ms1_list) - 1

                    return precursor_idx

                precursor_idx = search_precursor_scan(gpsm_list, ms1_list)

                def quant_one_gpsm(gpsm, ms1_list, precursor_idx, fq, need_mixed_check):
                    glymass = self.glycalc.calc_glycan_mass(gpsm.glycan)
                    modmass = self.pepcalc.calc_mod_mass_list(gpsm.seq, gpsm.mod)
                    bions = self.pepcalc.calc_b_ions(gpsm.seq, modmass)
                    pepmass = self.pepcalc.calc_pepmass_from_b(gpsm.seq, modmass, bions)
                    mono_mz = (
                        glymass + pepmass
                    ) / gpsm.charge + self.pepcalc.base_mass.mass_proton

                    # print(gpsm.spectrum, glymass, pepmass, mono_mz, gpsm.precursor_mz)
                    isotope_dist, mono_idx = self.isotope_calc.get_distribution(
                        gpsm.seq, gpsm.mod, gpsm.glycan, self.config.glyco_names
                    )
                    isotope_dist = np.array(isotope_dist)
                    isotope_mz = mono_mz * np.ones(len(isotope_dist)) + (
                        self.pepcalc.base_mass.mass_isotope / gpsm.charge
                    ) * (np.arange(len(isotope_dist)) - mono_idx)

                    if need_mixed_check:
                        Ac_mz = (
                            mono_mz
                            - np.arange(gpsm.glycan[fuc_idx] // 2, 0, -1)
                            * mixed_delta
                            / gpsm.charge
                        )

                    pre_idx = precursor_idx[gpsm.scan]
                    s = pre_idx
                    e = pre_idx
                    while s >= 0:
                        if not self._check_window(
                            ms1_list[s][-2],
                            ms1_list[s][-1],
                            gpsm.RT,
                            gpsm.scan,
                            self.RT_win,
                            self.Scan_win,
                        ):
                            break
                        s -= 1
                    while e < len(ms1_list):
                        if not self._check_window(
                            ms1_list[e][-2],
                            ms1_list[e][-1],
                            gpsm.RT,
                            gpsm.scan,
                            self.RT_win,
                            self.Scan_win,
                        ):
                            break
                        e += 1
                    if s < 0:
                        s = 0
                    if e >= len(ms1_list):
                        e = len(ms1_list) - 1

                    XIC_intens = []
                    RT_list = []
                    scan_list = []
                    XIC_Ac_intens = []
                    for peaklist, peakidx, RT, scan in ms1_list[s : e + 1]:
                        if need_mixed_check:
                            _, _, Ac_intens = spi.Match(
                                peaklist, peakidx, Ac_mz, self.tol, self.isppm
                            )
                        else:
                            Ac_intens = np.zeros(2)
                        idx, tol, intens = spi.Match(
                            peaklist, peakidx, isotope_mz, self.tol, self.isppm
                        )
                        # if scan == 19316:
                        # for mz,inten in peaklist:
                        # if int(mz) == 1112:
                        # print(peaklist[idx[0],:], intens[0], mz, mono_mz, inten)
                        XIC_intens.append(intens)
                        RT_list.append(RT)
                        scan_list.append(scan)
                        XIC_Ac_intens.append(Ac_intens)

                    XIC_intens = np.array(XIC_intens, dtype=float)
                    XIC_Ac_intens = np.array(XIC_Ac_intens, dtype=float)

                    used_iso = np.argmax(isotope_dist)
                    # used_iso = mono_idx
                    gpsm_idx = pre_idx - s

                    left, right, max_pos = self.quant.Retention_range(
                        XIC_intens[:, used_iso], gpsm_idx
                    )
                    # left, right = self.quant.Retention_range(XIC_intens[:,used_iso], gpsm_idx)

                    # left, right, max_pos = self.quant.side_trim(XIC_intens[:,used_iso], gpsm_idx, left, right, max_pos)

                    gpsm.mono_area = np.sum(XIC_intens[left : right + 1, mono_idx])
                    gpsm.peak_area = np.sum(XIC_intens[left : right + 1, :])

                    glycan = ""
                    glycan_composition = ""
                    if need_mixed_check:
                        ratio_Ac_intens = (
                            np.sum(XIC_Ac_intens[left : right + 1, :], axis=0)
                            / gpsm.mono_area
                        )
                        nonAcs = np.arange(ratio_Ac_intens.shape[0], 0, -1)[
                            ratio_Ac_intens < 0.1
                        ]
                        if len(nonAcs) > 0:
                            nAc = nonAcs[-1] - 1
                        else:
                            nAc = ratio_Ac_intens.shape[0]
                        correct_mono_area = 0
                        correct_peak_area = 0
                        if nAc != 0:
                            sum_Ac_intens = ratio_Ac_intens[-nAc:] * gpsm.mono_area
                            glycan = [gly for gly in gpsm.glycan]
                            glycan[fuc_idx] -= nAc * 2
                            glycan[ac_idx] += nAc
                            glycan_composition = "".join(
                                [
                                    "{}({})".format(name, gly) if gly != 0 else ""
                                    for gly, name in zip(
                                        glycan, self.config.glyco_names
                                    )
                                ]
                            )
                            glycan = " ".join([str(gly) for gly in glycan])
                            correct_mono_area = sum_Ac_intens[0]
                            correct_peak_area = gpsm.peak_area + np.sum(sum_Ac_intens)

                    # fqd.write("Spectrum=%s\nPeptide=%s\nModification=%s\nGlycan=%s\n"%(gpsm.spectrum, gpsm.seq, gpsm.mod, "".join(["%s(%d)"%(self.config.glyco_names[i],gpsm.glycan[i]) for i in range(len(gpsm.glycan)) if gpsm.glycan[i] > 0])))
                    # fqd.write("MS2_Scan=%d\n"%(gpsm.scan))
                    # fqd.write("MS2_RT=%.3f\n"%(gpsm.RT))
                    # fqd.write("Precursor_MZ=%.5f\n"%gpsm.precursor_mz)
                    # fqd.write("Precursor_Scan=%d\n"%ms1_list[pre_idx][-1])
                    # fqd.write("Precursor_RT=%.3f\n"%ms1_list[pre_idx][-2])
                    # fqd.write("Start_End_Position=%d\t%d\n"%(left, right))
                    # fqd.write("Start_End_RT=%.3f\t%.3f\n"%(RT_list[left], RT_list[right]))
                    # fqd.write("Start_End_Scan=%d\t%d\n"%(scan_list[left], scan_list[right]))
                    # fqd.write("Isotope_Num=%d\n"%len(isotope_dist))
                    # fqd.write("Theory_Isotope_MZ=%s\n"%("\t".join("%.5f"%isotope for isotope in isotope_mz)))
                    # fqd.write("Theory_Isotope_Intensity=%s\n"%("\t".join(["%.1f"%dist for dist in (isotope_dist/max(isotope_dist)*100)])))
                    # fqd.write("Extracted_MS1_RTs=\t%s\n"%("\t".join(["%.3f"%RT for RT in RT_list])))
                    # fqd.write("Extracted_MS1_Scans=\t%s\n"%("\t".join([str(scan) for scan in scan_list])))
                    # for i, intens in enumerate(XIC_intens.T):
                    # fqd.write("Intensities_Isotope#%d=\t%s\n"%(i, '\t'.join(['%.1f'%inten for inten in intens])))

                    fq.write("\t".join(gpsm.items))
                    # fq.write("\t%.5f"%(isotope_mz[mono_idx]))
                    fq.write("\t%g\t%g" % (gpsm.mono_area, gpsm.peak_area))
                    if not need_mixed_check:
                        fq.write("\t" + "\t".join(empty_correct))
                    else:
                        fq.write(
                            "\t%s\t%s\t%g\t%g"
                            % (
                                glycan,
                                glycan_composition,
                                correct_mono_area,
                                correct_peak_area,
                            )
                        )
                    fq.write("\n")

                for gpsm in gpsm_list:
                    need_check = False
                    if fuc_idx is not None and gpsm.glycan[fuc_idx] >= 2:
                        if ion_204_idx is not None:
                            inten_204 = float(gpsm.items[ion_204_idx])
                        else:
                            inten_204 = 0.001
                        marker_intens = []
                        if ion_274_idx is not None:
                            marker_intens.append(
                                float(gpsm.items[ion_274_idx]) / inten_204
                            )
                        if ion_292_idx is not None:
                            marker_intens.append(
                                float(gpsm.items[ion_292_idx]) / inten_204
                            )
                        marker_intens = np.array(marker_intens)

                        if np.all(marker_intens > 0.05) and np.any(marker_intens > 0.1):
                            need_check = True

                    quant_one_gpsm(gpsm, ms1_list, precursor_idx, fq, need_check)
                    cnt += 1
                    print(
                        "[XIC] Quantifying: %.1f%% (%d/%d)"
                        % (
                            float(cnt) / len(tmp_gpsm_list) * 100,
                            cnt,
                            len(tmp_gpsm_list),
                        ),
                        end="\r",
                    )

        fq.close()
        # fqd.close()
        for _reader in reader_dict.values():
            if _reader:
                _reader.close()
        # for _reader in pf2reader_dict.values():
        # if _reader: _reader.close()
        print("\n")


if __name__ == "__main__":
    import sys

    def ParseArgs(argv):
        arg_dict = {}
        i = 1
        while i < len(argv):
            arg_dict[argv[i]] = argv[i + 1]
            i += 2
        return arg_dict

    argd = ParseArgs(sys.argv)

    if len(sys.argv) < 3:
        print("Usage: XIC.exe -p pGlyco.cfg")
    else:
        cfg = argd["-p"]
        if "-w" in argd:
            window = int(argd["-w"])
        else:
            window = 21
        print("[XIC] Smoothing window = %d" % window)

        if "-m" in argd:
            method = argd["-m"]
        else:
            method = "savgol_filter"
        print("[XIC] Smoothing method = {}".format(method))

        xic = XICExtractor(cfg, method, window)

        if "-g" in argd:
            xic.RT_group = float(argd["-g"])
            print("[XIC] GPSMs are grouped by %.0f seconds" % xic.RT_group)
        if "-r" in argd:
            xic.PeakArea(pglyco_file=argd["-r"])
        else:
            xic.PeakArea(pglyco_file=None)
