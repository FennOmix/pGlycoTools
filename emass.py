from chemical import chem_dict
import numpy as np

class EMassSimple:
    def __init__(self, max_isotope_len):
        self.iso_dict = {}
        self.mono_dict = {}
        self.MAX_ISOTOPE_LEN = max_isotope_len
        for chem, isotope in chem_dict.items():
            tmp_dist = np.zeros(self.MAX_ISOTOPE_LEN)
            tmp_dist[0] = isotope[0][1]
            first_mass = isotope[0][0]
            mono_idx = 0
            max_inten = tmp_dist[0]
            for mass, inten in isotope[1:]:
                ith = round(mass - first_mass)
                if ith < self.MAX_ISOTOPE_LEN:
                    tmp_dist[ith] = inten
                    if max_inten < inten:
                        mono_idx = ith
                        max_inten = inten
            self.iso_dict[chem] = tmp_dist
            self.mono_dict[chem] = mono_idx

        self.zero_dist = np.zeros(self.MAX_ISOTOPE_LEN)
        self.zero_dist[0] = 1

    def one_chem_dist(self, chem, n):
        if n == 0:
            return self.zero_dist, 0
        elif n == 1:
            return self.iso_dict[chem], self.mono_dict[chem]
        tmp_dist, mono_idx = self.one_chem_dist(chem, int(n / 2))
        tmp_dist, mono_idx = self.convolution(tmp_dist, mono_idx, tmp_dist, mono_idx)
        if n % 2 == 1:
            return self.convolution(
                tmp_dist, mono_idx, self.iso_dict[chem], self.mono_dict[chem]
            )
        else:
            return tmp_dist, mono_idx

    def convolution(self, d1, mono_idx1, d2, mono_idx2):
        mono_idx = mono_idx1 + mono_idx2
        ret = np.zeros(self.MAX_ISOTOPE_LEN * 2 - 1)
        for i in range(len(d1)):
            for j in range(len(d2)):
                ret[i + j] += d1[i] * d2[j]

        # keep top-k peaks around mono
        mono_left = mono_idx - 1
        mono_right = mono_idx + 1
        while (
            mono_left >= 0
            and mono_right < len(ret)
            and (mono_right - mono_left - 1) < self.MAX_ISOTOPE_LEN
        ):
            if ret[mono_right] >= ret[mono_left]:
                mono_right += 1
            else:
                mono_left -= 1
        while mono_left == -1 and (mono_right - mono_left - 1) < self.MAX_ISOTOPE_LEN:
            mono_right = self.MAX_ISOTOPE_LEN
        while (
            mono_right == len(ret)
            and (mono_right - mono_left - 1) < self.MAX_ISOTOPE_LEN
        ):
            mono_left = (
                self.MAX_ISOTOPE_LEN - 2
            )  # self.MAX_ISOTOPE_LEN*2-1-self.MAX_ISOTOPE_LEN - 1

        return ret[mono_left + 1 : mono_right], mono_idx - mono_left - 1

    # list formula, for example: [('H', 2), ('O', 1)]
    def get_distribution(self, formula):
        ret = self.zero_dist
        ret_mono = 0
        for chem, n in formula:
            _one_chem_dist, one_chem_mono = self.one_chem_dist(chem, n)
            ret, ret_mono = self.convolution(
                ret, ret_mono, _one_chem_dist, one_chem_mono
            )
        return ret, ret_mono
